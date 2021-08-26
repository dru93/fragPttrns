library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicAlignments)
library(GenomicRanges)
library(dplyr)
library(readr)
library(httr)

# get fragment ratio of one sample in bins of AB GRange object
getFragmentRatioFromBAM = function(sampleFilename, AB,
                                   bins5Mb = F, print = T, trimFRs = F, mapqFilt = 30, FRthres = 150){
  # sampleFilename: full path name of BAM file, AB: empty 100kb bins GRange object, getHg19bins100kb
  # bins5Mb: aggregate to 5Mb windows, trimFRs: trim from to 
  # mapqFilt: MAPQ filter of BAM file, FRthres: short/long fragment desicion threshold
  gc.correct = function(coverage, bias){ # correct for GC content
    i = seq(min(bias, na.rm = TRUE), max(bias, na.rm = TRUE), by = 0.001)
    coverage.trend = loess(coverage ~ bias)
    coverage.model = loess(predict(coverage.trend, i) ~ i)
    coverage.pred = predict(coverage.model, bias)
    coverage.corrected = coverage - coverage.pred + median(coverage, na.rm = T)
  }
  Mode = function(x){ # get mode of vector
    ux = unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  # retrieve files basename
  sampleName = basename(sampleFilename)
  # set bam parameters to keep and read alignment pairs of BAM file
  if (print) cat(paste0("sample:", sampleName, "\nReading BAM file\n\n"))
  param = ScanBamParam(mapqFilter = mapqFilt,
                       flag = scanBamFlag(isDuplicate = F, isSecondaryAlignment = F, isUnmappedQuery = F))
  alignments = readGAlignmentPairs(sampleFilename, param = param, index = sampleFilename)
  fragments = granges(keepSeqlevels(alignments, paste0("chr", 1:22), pruning.mode = "coarse"), on.discordant.seqnames = "drop")
  rm(alignments, param) # help RAM
  gc()
  # filter outliers
  w.all = width(fragments)
  q.all = quantile(w.all, c(0.001, 0.999))
  fragments = fragments[which(w.all > q.all[1] & w.all < q.all[2])]
  unstrfragments = unstrand(fragments) # remove strand information
  # compute gc content
  if (print) cat(paste0("sample:", sampleName, "\nComputing GC content\n\n"))
  gcs = biovizBase::GCcontent(Hsapiens, unstrfragments)
  fragments$gc = gcs
  # filter DUKE blacklisted regions
  filters.hg19 = getHg19blacklistedRegions() 
  fragments = fragments[-queryHits(findOverlaps(fragments, filters.hg19))]
  w.all = width(fragments) # calculate fragment size
  # keep only fragments with size between 100-220
  if (trimFRs){
    fragments = fragments[which(w.all >= 100 & w.all <= 220)]
    w = width(fragments)
  } else {
    w = w.all
  }
  frag.list = split(fragments, w) # list of fragments by fragment size
  # get fragment size in AB bins
  if (print) cat(paste0("sample:", sampleName, "\nRetrieving frament size in bins\n\n"))
  # counts is a count matrix with y: unique fragment sizes as columns and x: AB bins
  counts = sapply(frag.list, function(x) countOverlaps(AB, x))
  if(min(w) > 100) {
    m0 = matrix(0, ncol = min(w) - 100, nrow = nrow(counts), dimnames = list(rownames(counts), 100:(min(w)-1)))
    counts = cbind(m0, counts)
  }
  # get gc content of 100kb bins
  if (print) cat(paste0("sample:", sampleName, "\nRetrieving bins GC content\n\n"))
  olaps = findOverlaps(fragments, AB)
  bin.list = split(fragments[queryHits(olaps)], subjectHits(olaps))
  bingc = rep(NA, length(bin.list))
  bingc[unique(subjectHits(olaps))] = sapply(bin.list, function(x) mean(x$gc))
  # get ratio of shot long fragment size
  colIdxThres = which(colnames(counts) == FRthres)
  short = rowSums(counts[,1:colIdxThres]) 
  long = rowSums(counts[,(colIdxThres + 1):ncol(counts)])
  ratio = short/long
  # replace inf with NA @ ratio in case of division with 0 (long = 0)
  ratio = replace(ratio, is.infinite(ratio), NA)
  # agregate eveything to AB and correct for GC content
  if (print) cat(paste0("sample:", sampleName, "\nCorrecting for GC content\n\n"))
  AB$short = short
  AB$long = long
  AB$short.corrected = gc.correct(short, bingc)
  AB$long.corrected = gc.correct(long, bingc)
  AB$nfrags.corrected = gc.correct(short + long, bingc)
  AB$ratio.corrected = gc.correct(ratio, bingc)
  AB$mode = Mode(w)
  AB$mean = round(mean(w), 2)
  AB$median = median(w)
  AB$quantile.25 = quantile(w, 0.25)
  AB$quantile.75 = quantile(w, 0.75)
  AB$frag.gc = bingc
  # add fragment size for each of the columns/fragment lengths
  for(i in 1:ncol(counts)) elementMetadata(AB)[,colnames(counts)[i]] = counts[,i]
  df.fr = as.data.frame(AB)
  df.fr$arm = as.factor(df.fr$arm)
  if (!bins5Mb) {
    return(df.fr)
  } else {
    # combine variable annotates rows to be combined 5Mb = (50*100kb)
    df.fr2 = df.fr %>% group_by(arm) %>%
      mutate(combine = ifelse(grepl("p", arm), ceiling((1:length(arm))/50),
                              ceiling(rev((1:length(arm))/50) )))
    # merging to 5Mb bins
    df.fr3 = df.fr2 %>% group_by(seqnames, arm, combine) %>%
      summarize(short = mean(short),
                long = mean(long),
                short.corrected2 = mean(short.corrected),
                long.corrected2 = mean(long.corrected),
                eigen = mean(eigen),
                gc = mean(frag.gc),
                ratio.corrected = mean(ratio.corrected),
                nfrags.corrected2 = mean(nfrags.corrected),
                domain = names(which.max(table(domain))),
                short.var = var(short.corrected, na.rm = T),
                long.var = var(long.corrected),
                nfrags.var = var(nfrags.corrected),
                mode_size = unique(mode),
                mean_size = unique(mean),
                median_size = unique(median),
                q25_size = unique(quantile.25),
                q75_size = unique(quantile.75),
                start = start[1],
                end = rev(end)[1],
                binsize = n(), .groups = "drop_last")
    # keep 5Mb bins (50*100kb)
    df.fr3 = df.fr3 %>% filter(binsize == 50)
    # add bin and ratio.centered columns, remove combine column, correct names
    df.fr3$bin = c(1:nrow(df.fr3))
    df.fr3$ratio.centered = scale(df.fr3$ratio.corrected, center = T, scale = T)
    df.fr3$combine = NULL
    names(df.fr3)[which(names(df.fr3) == "short.corrected2")] = "short.corrected"
    names(df.fr3)[which(names(df.fr3) == "long.corrected2")] = "long.corrected"
    names(df.fr3)[which(names(df.fr3) == "nfrags.corrected2")] = "nfrags.corrected"
    return(df.fr3)
  }
}

# get hg19 gaps
getHg19gaps = function(){
  genome = "hg19"
  mySession = browserSession()
  genome(mySession) = genome
  # get hg19 gaps
  gaps = getTable(ucscTableQuery(mySession, track = "gap"))
  # convert gaps to GRanges object and keep 1-22 & X, Y chromosomes
  gaps.hg19 = GRanges(gaps$chrom, IRanges(gaps$chromStart, gaps$chromEnd), type = gaps$type)
  gaps.hg19 = keepSeqlevels(gaps.hg19, paste0("chr", c(1:22, "X", "Y")), pruning.mode = "coarse")
  seqinfo(gaps.hg19) = seqinfo(Hsapiens)[seqlevels(gaps.hg19),]
  return(gaps.hg19)
}

# get Duke blacklisted regions
getHg19blacklistedRegions = function(){
  url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"
  blacklisted.file = httr::content(GET(url))
  blacklisted.tib = read_tsv(gzcon(rawConnection(blacklisted.file)),
                             col_names = c("seqnames", "start", "end", "name", "score", "."), col_types = cols())
  blacklisted.tib = blacklisted.tib %>% mutate(start = start+1)
  filters.hg19 = makeGRangesFromDataFrame(blacklisted.tib, keep.extra.columns = TRUE)
  filters.hg19 = keepSeqlevels(filters.hg19, paste0("chr", c(1:22, "X", "Y")), pruning.mode = "coarse")
  seqinfo(filters.hg19) = seqinfo(Hsapiens)[seqlevels(filters.hg19),]
  return(filters.hg19)
}

# load 100kb empty bins
getHg19bins100kb = function(){
  # load empty 100kb bin genome fragments
  ABurl = RCurl::getURL('https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt', ssl.verifyhost = FALSE, ssl.verifypeer = FALSE)
  AB = read.table(textConnection(ABurl), header = TRUE)
  AB = makeGRangesFromDataFrame(AB, keep.extra.columns = TRUE)
  # filter AB
  chromosomes = GRanges(paste0("chr", 1:22), IRanges(0, seqlengths(Hsapiens)[1:22]))
  gaps.hg19 = getHg19gaps()
  tcmeres = gaps.hg19[grepl("centromere|telomere", gaps.hg19$type)]
  arms = GenomicRanges::setdiff(chromosomes, tcmeres)
  arms = arms[-c(25,27,29,41,43)]
  armlevels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
                "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
                "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
                "19p", "19q","20p","20q","21q","22q")
  arms$arm = armlevels
  AB = AB[-queryHits(findOverlaps(AB, gaps.hg19))]
  AB = AB[queryHits(findOverlaps(AB, arms))]
  AB$arm = armlevels[subjectHits(findOverlaps(AB, arms))]
  seqinfo(AB) = seqinfo(Hsapiens)[seqlevels(seqinfo(AB))]
  AB = trim(AB)
  # get AB gc content
  AB$gc = biovizBase::GCcontent(Hsapiens, AB)
  ## These bins had no coverage
  AB = AB[-c(8780, 13665)]
  return(AB)
}