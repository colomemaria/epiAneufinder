#' Generate count matrix
#' 
#' @param reads BamFileList, Grange object or dgCMatrix object dependent 
#'              on the original input file type (see main wrapper function)
#' @param windows Binned genome
#' @param by barcode information
#' @param minFrags Minimum number of fragments a barcode must contain to be counted as a cell
#' @param mapqFilter Filter bam files after a certain mapq value
#' @export
generateCountMatrix <- function(reads, windows, by=NULL, minFrags=NULL, mapqFilter=NULL){

  #Keep only regions in filtered chromosomes
  windows   <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")

  #Count Insertions in windows
  message("Getting Counts...")
  counts <- countInsertions(reads, windows, by="barcode", 
                            minFrags=minFrags,mapqFilter=mapqFilter)

  #Keep only regions with less than 0.1% N
  keep <- which(windows$N < 0.001)
  windowSummary <- windows[keep,]
  countSummary <- counts[keep,]

  se <- SummarizedExperiment(assays = list(counts = as.matrix(countSummary)),
                             rowRanges = windowSummary)
  colnames(se) <- colnames(counts)

  return(se)
}
