#' @import GenomicFeatures
#' @export

addExpressionFactor <- function(bins, gene.annotation=NULL) {
  UseMethod("addExpressionFactor")
}

# Looking at genes
addExpressionFactor.GRanges <- function(bins, gene.annotation=NULL) {
  txdb <- getFromNamespace(gene.annotation, ns=gene.annotation)
  seqlevelsStyle(txdb) <- seqlevelsStyle(bins)[1]
  genes <- sort(keepStandardChromosomes(genes(txdb), pruning.mode = 'coarse'))
  bins$genecoverage <- countOverlaps(bins,genes)
  # bins
}

addExpressionFactor.list <- function(bins, gene.annotation=NULL) {
  lapply(bins, addExpressionFactor, txdb)
}

addExpressionFactor.default <- function(bins, gene.annotation=NULL) {
  stop("Do not know how to get expression factor for type ", sQuote(class(bins)))
}
