#' @export
#' @param reads bamfile or fragments file of sequencing information
#' @param windows Binned genome

generateCountMatrix <- function(reads, windows, by=NULL, minFrags=NULL){

  #Keep only regions in filtered chromosomes
  windows   <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")

  #Count Insertions in windows
  message("Getting Counts...")
  counts <- countInsertions(reads, windows, by="barcode", minFrags=minFrags)

  #Keep only regions with less than 0.1% N
  keep <- which(windows$N < 0.001)
  windowSummary <- windows[keep,]
  countSummary <- counts[keep,]

  se <- SummarizedExperiment(assays = list(counts = as.matrix(countSummary)), rowRanges = windowSummary)
  colnames(se) <- colnames(counts)

  return(se)
}
