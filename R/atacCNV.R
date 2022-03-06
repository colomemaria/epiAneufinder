#' Wrapper function for the \code{\link{atacCNV}} package
#'
#' @param input Folder with bam files or a fragments.tsv file
#' @param outdir Path to output directory
#' @param blacklist Bed file with blacklisted regions
#' @param windowSize Size of the window (Reccomended for sparse data - 1e6)
#' @param genome String containing name of BS.genome object. Necessary for GC correction. Default: "BSgenome.Hsapiens.UCSC.hg38"
#' @param test One of AD or KS. Anderson-Darling or Kolmogorov-Smirnov. Default: "AD"
#' @param reuse.existing Logical. False removes all the files in the outdir and recomputes everything.
#' @param exclude String of chromosomes to exclude. Example: c('chrX','chrY','chrM')
#' @param readout String. One of "ATAC", "BS". Default: "ATAC"
#' @param uq Upper quantile. Default: 0.1
#' @param lq Lower quantile. Default: 0.9
#' @param title_karyo String. Title of the output karyogram
#' @param minFrags Integer. Minimum number of reads for a cell to pass. Only required for fragments.tsv file. Default: 20000
#' @param gene.annotation String. Txdb object for gene presence correction. Set to NULL to not correct. Default: NULL
#' @param threshold_blacklist_bins Blacklist a bin if more than the given ratio of cells have zero reads in the bin. Default: 0.85
#' @param ncores Number of cores for parallelization. Default: 4
#' @param minsize Integer. Resolution at the level of ins. Default: 1. Setting it to higher numbers runs the algorithm faster at the cost of resolution
#' @param k Integer. Find 2^k segments per chromosome
#' @import stats
#' @import GenomicRanges
#' @import plyranges
#' @import S4Vectors
#' @import BiocGenerics
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import data.table
#' @import Rsamtools
#' @import GenomeInfoDb
#' @import Biostrings
#' @import Matrix
#' @import ggplot2
#' @import cowplot
#' @import SummarizedExperiment
#' @import R.utils
#' @import magrittr
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom parallel mclapply
#' @return \code{NULL}
#' @author Akshaya Ramakrishnan
#' @export


# Imports:
# GenomicAlignments (>= 1.22.1),
# Matrix (>= 1.2.18),
# SummarizedExperiment (>= 1.16.1),
# data.table (>= 1.12.8),
# plyranges (>= 1.6.10),
# Rsamtools (>= 2.2.1),
# GenomeInfoDb (>= 1.22.0),
# BSgenome.Hsapiens.UCSC.hg38 (>= 1.4.1),
# GenomicRanges (>= 1.38.0),
# stats (>= 3.6.1),
# Biostrings (>= 2.54.0),
# magrittr (>= 1.5),
# BiocGenerics (>= 0.32.0),
# S4Vectors (>= 0.24.2)

atacCNV <- function(input, outdir, blacklist, windowSize, genome="BSgenome.Hsapiens.UCSC.hg38",
                    test='AD', reuse.existing=FALSE, exclude=NULL, readout="ATAC",
                    uq=0.9, lq=0.5, somyl=0.1, somyu=0.8, title_karyo=NULL, minFrags = 20000,
                    gene.annotation=NULL, threshold_blacklist_bins=0.85, ncores=4, minsize=1, k=3){

  if(reuse.existing==FALSE){
    print("Removing old file from the output folder")
    file.remove(list.files(outdir, full.names=TRUE))
  }

  if(!file.exists(file.path(outdir,"count_summary.rds"))) {
    blacklist <- read_bed(blacklist)
    windows <- makeWindows(genome = genome, blacklist = blacklist, windowSize, exclude=exclude)

    if(file_test("-d", input)){
      print("Obtaining bam file list")
      bamfiles <- Rsamtools::BamFileList(list.files(input, pattern = ".bam$", full.names = TRUE), yieldSize=100000)
      # print(bamfiles)
      counts <- generateCountMatrix(bamfiles, windows)
    }
    else if(file_test("-f", input)){
      if(grepl("\\.tsv$|\\.tsv.gz$", input)){
        print("Obtaining the fragments tsv file")
        file_fragments <- fread(input)
        colnames(file_fragments) <- c('seqnames','start','end','barcode','pcr')
        fragments <- as_granges(file_fragments)
        print(head(fragments))
        print(names(fragments))
      } else if(grepl("\\.bed$", input)){
        print("Obtaining the fragments bed file")
        fragments <- read_bed(input)
        names(mcols(fragments)) <- 'barcode'
      } else{
        stop("Please provide a fragments .tsv/.bed or a path to the directory containing all the bam files")
      }
      counts <- generateCountMatrix(fragments, windows, by="barcode", minFrags = minFrags)
    }
    saveRDS(counts, file.path(outdir,"count_summary.rds"))
  }

  counts <- readRDS(file.path(outdir,"count_summary.rds"))
  peaks <- as.data.table(assays(counts)$counts)
  colnames(peaks) <- paste0('cell-', colnames(peaks))
  rowinfo <- as.data.table(rowRanges(counts))
  peaks <- cbind(rowinfo, peaks)

  if(!file.exists(file.path(outdir,"counts_gc_corrected.rds"))) {
    message("Correcting for GC bias...")
    corrected_counts <- peaks[, mclapply(.SD, function(x) {
      fit <- stats::loess(x ~ peaks$GC)
      correction <- mean(x) / fit$fitted
      as.integer(round(x * correction))
    }, mc.cores = ncores), .SDcols = patterns("cell-")]
    saveRDS(corrected_counts, file.path(outdir,"counts_gc_corrected.rds"))
  }

  corrected_counts <- readRDS(file.path(outdir,"counts_gc_corrected.rds"))
  peaks <- cbind(rowinfo, corrected_counts)

  zeroes_per_bin <- peaks[, rowSums(.SD==0), .SDcols = patterns("cell-")]
  ncells <- length(grep("cell-", colnames(peaks)))
  peaks <- peaks[zeroes_per_bin<(threshold_blacklist_bins*ncells)]

  if(!file.exists(file.path(outdir,"results_gc_corrected.rds"))) {
    clusters_ad <- peaks[, mclapply(.SD, function(x) {
      peaksperchrom <- split(x, peaks$seqnames)
      print("Calculating distance AD")
      results <- lapply(peaksperchrom, function(x2) {
        getbp(x2, k = k, minsize = minsize, test=test)
      })
    }, mc.cores = ncores), .SDcols = patterns("cell-")]
    saveRDS(clusters_ad, file.path(outdir, "results_gc_corrected.rds"))
  }
  print("Successfully identified breakpoints")

  if(!file.exists(file.path(outdir,"cnv_calls.rds"))) {
    names_seq <- levels(peaks$seqnames)

    clusters_ad <- readRDS(file.path(outdir,"results_gc_corrected.rds"))
    breakpoints <- lapply(clusters_ad, function(x) { lapply(x,'[[', 1) })
    distances <- lapply(clusters_ad, function(x) { lapply(x,'[[', 2) })
    result.dt <- Map(function(bp, dist){
      names(bp) <- names_seq
      names(dist) <- names_seq
      dtlist <- Map(function(per_chr_bp, per_chr_dist){
        data.table(per_chr_bp, per_chr_dist)
      }, bp, dist)
      per_cell_dt <- rbindlist(dtlist, idcol = 'Chr')
      per_cell_dt <- per_cell_dt[order(-per_chr_dist)]
    }, breakpoints, distances)

    pruned_result.dt <- lapply(result.dt, function(x){
      threshold_dist_values(x)
    })
    print("Successfully discarded irrelevant breakpoints")


    clusters_pruned <- as.data.table(Map(function(seq_data, bp){
      per_chrom_seq_data <- split(seq_data, peaks$seqnames)
      per_chrom_bp <- split(bp, bp$Chr)
      clusters_per_chrom <- Map(function(seq_data2, bp2){
        if(is.null(bp2)){
          return(rep(1, length(seq_data2)))
        } else {
          bp_to_cluster <- sort(c(1,length(seq_data2)+1,bp2$per_chr_bp))
          return(rep(1:length(diff(bp_to_cluster)),diff(bp_to_cluster)))
        }
      }, per_chrom_seq_data[names_seq], per_chrom_bp[names_seq])
      cl <- 0
      clusters <- vector()
      for(item in clusters_per_chrom){
        item <- item + cl
        clusters <- append(clusters, item)
        cl <- cl + length(unique(item))
      }
      return(clusters)
    },  peaks[, .SD, .SDcols = patterns("cell-")], pruned_result.dt))

    if(readout=="ATAC"){
      somies_ad <- Map(function(seq_data,cluster) {
        assign_gainloss(seq_data, cluster, uq=uq, lq=lq)
      }, peaks[, .SD, .SDcols = patterns("cell-")], clusters_pruned)
      print("Successfully assigned gain-loss")
    }
    if(readout=="BS"){
      somies_ad <- Map(function(seq_data,cluster) {
        assign_somy(seq_data, cluster, uq=uq, lq=lq, somyl=somyl, somyu=somyu)
      }, peaks[, .SD, .SDcols = patterns("cell-")], clusters_pruned)
      print("Successfully assigned gain-loss")
    }
    saveRDS(somies_ad, file.path(outdir, "cnv_calls.rds"))
  }

  somies_ad <- readRDS(file.path(outdir,"cnv_calls.rds"))
  if(is.null(title_karyo)){
    title_karyo <- basename(outdir)
  }
  if(readout=="ATAC"){
    plot_karyo_gainloss(somies_ad = somies_ad, outdir = outdir, peaks = peaks, uq, lq, title_karyo)
  }
  if(readout=="BS"){
    plot_karyo(somies_ad = somies_ad, outdir = outdir, peaks = peaks, uq, lq, somyl, somyu, title_karyo)
  }
  print("Successfully plotted karyogram")
}
