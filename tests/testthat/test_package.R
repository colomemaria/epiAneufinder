input <- "/Users/akshaya/Work/atacCNV/data/ATAC/HCT116"
outdir <- "/Users/akshaya/Work/atacCNV/workspace/ad_bam/package_test/"
blacklist <- "/Users/akshaya/Work/hg38.blacklist.bed"
windowSize <- 1e6
atacCNV(input, outdir, blacklist, windowSize)

