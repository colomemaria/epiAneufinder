#' @export
#' @param genome String containing name of BS.genome object. Necessary for GC correction. Default: "BSgenome.Hsapiens.UCSC.hg38"
#' @param blacklist Granges object with information about blacklisted regions in the genome.
#' @param windowSize Size in basepairs for a bin
#' @param exclude Chromosomes to exclude

makeWindows <- function(genome, blacklist, windowSize, exclude = NULL){
  genome <- getFromNamespace(genome, ns=genome)
  windows <- tileGenome(seqlengths = seqlengths(genome), tilewidth = windowSize, cut.last.tile.in.chrom = TRUE)
  windows <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
  windows <- dropSeqlevels(windows, exclude, pruning.mode = 'coarse')
  mcols(windows)$wSeq <- as.character(seqnames(windows))
  mcols(windows)$wStart <- BiocGenerics::start(windows)
  mcols(windows)$wEnd <- BiocGenerics::end(windows)
  message("Subtracting Blacklist...")
  overlaps <- findOverlaps(windows, blacklist)
  idx <- setdiff(1:length(windows), S4Vectors::queryHits(overlaps))
  windowsBL <- windows[idx]
  names(windowsBL) <- paste0("w",seq_along(windowsBL))
  mcols(windowsBL)$name <- names(windowsBL)
  message("Adding Nucleotide Information...")
  windowSplit <- split(windowsBL, as.character(seqnames(windowsBL)))
  windowNuc <- lapply(seq_along(windowSplit), function(x){
    message(sprintf("%s of %s", x, length(windowSplit)))
    chrSeq <- Biostrings::getSeq(genome,names(windowSplit)[x])
    grx <- windowSplit[[x]]
    aFreq <- Biostrings::alphabetFrequency(Biostrings::Views(chrSeq, ranges(grx)))
    mcols(grx)$GC <- rowSums(aFreq[, c("G","C"),drop=FALSE]) / rowSums(aFreq)
    mcols(grx)$AT <- rowSums(aFreq[, c("A","T"),drop=FALSE]) / rowSums(aFreq)
    return(grx)
  }) %>% GRangesList %>% unlist %>% sortSeqlevels %>% sort
  windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
  print("Finished making windows successfully")
  windowNuc
}
