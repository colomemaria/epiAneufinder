#' Wrapper function for the \code{\link{epiAneufinder}} package
#'
#' @param input Folder with bam files, a fragments.tsv/bed file or a folder with
#'              a count matrix (required files: matrix.mtx(.gz), barcodes.tsv(.gz) and peaks.bed(.gz))
#' @param outdir Path to output directory
#' @param blacklist Bed file with blacklisted regions
#' @param windowSize Size of the window (Reccomended for sparse data - 1e6)
#' @param genome String containing name of BS.genome object. Necessary for GC correction. Default: "BSgenome.Hsapiens.UCSC.hg38"
#' @param test One of AD or KS. Anderson-Darling or Kolmogorov-Smirnov. Default: "AD"
#' @param reuse.existing Logical. False removes all the files in the outdir and recomputes everything.
#' @param exclude String of chromosomes to exclude. Example: c('chrX','chrY','chrM')
#' @param uq Upper quantile. Default: 0.1
#' @param lq Lower quantile. Default: 0.9
#' @param title_karyo String. Title of the output karyogram
#' @param minFrags Integer. Minimum number of reads for a cell to pass. Only required for fragments.tsv file. Default: 20000
#' @param mapqFilter Filter bam files after a certain mapq value
#' @param threshold_cells_nbins Keep only cells that have more than a certain percentage of non-zero bins
#' @param selected_cells Additional option for filtering the input, either NULL or a file with barcodes of cells to keep (one barcode per line, no header)
#' @param threshold_blacklist_bins Blacklist a bin if more than the given ratio of cells have zero reads in the bin. Default: 0.85
#' @param ncores Number of cores for parallelization. Default: 4
#' @param minsize Integer. Resolution at the level of ins. Default: 1. Setting it to higher numbers runs the algorithm faster at the cost of resolution
#' @param k Integer. Find 2^k segments per chromosome
#' @param minsizeCNV Integer. Number of consecutive bins to constitute a possible CNV
#' @param save_removed_regions Boolean variable. Option to save regions that were filtered out in a file removed_regions.tsv.
#' @param gc_correction Boolean variable. Option to perform GC correction. Strongly advised to set to TRUE unless the input matrix is already GC corrected!!!
#' @param plotKaryo Boolean variable. Whether the final karyogram is plotted at the end
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

epiAneufinder <- function(input, outdir, blacklist, windowSize, genome="BSgenome.Hsapiens.UCSC.hg38",
                    test='AD', reuse.existing=FALSE, exclude=NULL,
                    uq=0.9, lq=0.1, title_karyo=NULL, minFrags = 20000, mapqFilter=10,
                    threshold_cells_nbins=0.05,selected_cells=NULL,
                    threshold_blacklist_bins=0.85,
                    ncores=4, minsize=1, k=4, 
                    minsizeCNV=0,save_removed_regions=FALSE,gc_correction=TRUE,
                    plotKaryo=TRUE){

  outdir <- file.path(outdir, "epiAneufinder_results")
  dir.create(outdir,recursive=TRUE)
  
  if(reuse.existing==FALSE){
    message("Removing old file from the output folder")
    file.remove(list.files(outdir, full.names=TRUE))
  }

  #Restrict the number of cores inside BLAS and data table
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  data.table::setDTthreads(1)
  
  # ----------------------------------------------------------------------------
  # Creating the count matrix
  # ----------------------------------------------------------------------------
  
  if(!file.exists(file.path(outdir,"count_summary.rds"))) {
    blacklist <- read_bed(blacklist)
    windows <- makeWindows(genome = genome, blacklist = blacklist, windowSize, exclude=exclude)

    if(file_test("-d", input)){
      if(any(file.exists(file.path(input,c("matrix.mtx","matrix.mtx.gz"))))){
        message("Obtaining count matrix and reformting it to specified windows")
        peak_matrix<-readCountMatrix(input)
        
        counts <- generateCountMatrix(peak_matrix, windows)
      } else {
        message("Obtaining bam file list")
        bamfiles <- Rsamtools::BamFileList(list.files(input, pattern = ".bam$", full.names = TRUE), yieldSize=100000)
        # print(bamfiles)
        counts <- generateCountMatrix(bamfiles, windows,mapqFilter=mapqFilter)
      }
    }
    else if(file_test("-f", input)){
      if(grepl("\\.tsv$|\\.tsv.gz$", input)){
        message("Obtaining the fragments tsv file")
        file_fragments <- fread(input)
        colnames(file_fragments) <- c('seqnames','start','end','barcode','pcr')
        fragments <- as_granges(file_fragments)
        #print(head(fragments))
      } else if(grepl("\\.bed$", input)){
        message("Obtaining the fragments bed file")
        fragments <- read_bed(input)
        names(mcols(fragments)) <- 'barcode'
      } else{
        stop("Please provide a fragments .tsv/.bed or a path to the directory containing all the bam files")
      }
      counts <- generateCountMatrix(fragments, windows, by="barcode", minFrags = minFrags)
    }
    
    #Check that the count matrix really has entries
    if(ncol(counts)==0){
      stop("The created count matrix is empty. Please check input files and filtering options.")
    }
    
    message(paste("Count matrix with",ncol(counts),"cells and",nrow(counts),"windows",
                "has been generated and will be saved as count_summary.rds"))
    saveRDS(counts, file.path(outdir,"count_summary.rds"))
  }

  counts <- readRDS(file.path(outdir,"count_summary.rds"))
  peaks <- as.data.table(assays(counts)$counts)
  
  # ----------------------------------------------------------------------------
  # Filtering cells 1) based on barcode file (if provided) and
  # 2) based on too many zero windows and
  # 3) filter windows without enough coverage
  # ----------------------------------------------------------------------------
  
  # Filter cells based on a barcode file if provided
  if(! is.null(selected_cells)){
    cells_select<-read.table(selected_cells)
    
    #Filter for cells which are in the count matrix
    cells_select<-cells_select[cells_select$V1 %in% colnames(peaks),,drop=FALSE]
    
    peaks <- peaks[,cells_select$V1,with=FALSE]
    
    message(paste("Filtering cell based on additionally provided barcode file,",
                ncol(peaks),"cells remain."))
  }

  # Exclude cells that have no signal in most bins
  zeroes_per_cell<-colSums(peaks==0)
  peaks <- peaks[,zeroes_per_cell<((1-threshold_cells_nbins)*nrow(peaks)),with=FALSE]
  message(paste("Filtering cell without enough covered windows,",
              ncol(peaks),"cells remain."))
  
  # Add row information (window location)
  colnames(peaks) <- paste0('cell-', colnames(peaks))
  rowinfo <- as.data.table(rowRanges(counts))
  peaks <- cbind(rowinfo, peaks)
  
  #Exclude bins (=windows) with too little signal
  zeroes_per_bin <- peaks[, rowSums(.SD==0), .SDcols = patterns("cell-")]
  ncells <- length(grep("cell-", colnames(peaks)))
  
  # Save which bins will be removed in the next step as a separate file
  if(save_removed_regions){
    write.table(file=file.path(outdir,"removed_regions.tsv"),
                peaks[zeroes_per_bin>=(threshold_blacklist_bins*ncells),
                      c("seqnames","start","end")],
                sep="\t",quote=FALSE,row.names=FALSE)
  }

  # Exclude bins that have no signal in most cells
  peaks <- peaks[zeroes_per_bin<(threshold_blacklist_bins*ncells)]
  rowinfo<-rowinfo[zeroes_per_bin<(threshold_blacklist_bins*ncells)]
  
  #Drop factor levels of empty chromosomes
  peaks$seqnames<-droplevels(peaks$seqnames)
  
  message(paste("Filtering empty windows,",nrow(peaks),"windows remain."))
  
  # ----------------------------------------------------------------------------
  # GC correction
  # ----------------------------------------------------------------------------
    
  if(!file.exists(file.path(outdir,"counts_gc_corrected.rds"))) {
    
    if(gc_correction){
      
      message("Correcting for GC bias using a LOESS fit ...")
      
      corrected_counts <- peaks[, mclapply(.SD, function(x) {
        # LOESS correction for GC
        fit <- stats::loess(x ~ peaks$GC)
        correction <- mean(x) / fit$fitted
        as.integer(round(x * correction))
      }, mc.cores = ncores), .SDcols = patterns("cell-")]
    } else {
      
      message("Warning: skipping the GC correction!")
      
      corrected_counts <- peaks[,.SD,.SDcols = patterns("cell-")]
    }

    saveRDS(corrected_counts, file.path(outdir,"counts_gc_corrected.rds"))
  
  } 
  
  # ----------------------------------------------------------------------------
  # Estimating breakpoints
  # ----------------------------------------------------------------------------
  
  corrected_counts <- readRDS(file.path(outdir,"counts_gc_corrected.rds"))
  peaks <- cbind(rowinfo, corrected_counts)
  
  if(!file.exists(file.path(outdir,"results_gc_corrected.rds"))) {
    
    message("Calculating distance AD")
    
    clusters_ad <- peaks[, mclapply(.SD, function(x) {
      peaksperchrom <- split(x, peaks$seqnames)
      results <- lapply(peaksperchrom, function(x2) {
        getbp(x2, k = k, minsize = minsize, test=test,minsizeCNV=minsizeCNV)
      })
      return(results)
    }, mc.cores = ncores), .SDcols = patterns("cell-")]
    saveRDS(clusters_ad, file.path(outdir, "results_gc_corrected.rds"))
  }
  message("Successfully identified breakpoints")

  # ----------------------------------------------------------------------------
  # Pruning break points and annotating CNV status of each segment
  # ----------------------------------------------------------------------------
  
  if(!file.exists(file.path(outdir,"cnv_calls.rds"))) {
    names_seq <- levels(peaks$seqnames)

    clusters_ad <- readRDS(file.path(outdir,"results_gc_corrected.rds"))
    breakpoints <- lapply(clusters_ad, function(x) { lapply(x,'[[', 1) })
    distances <- lapply(clusters_ad, function(x) { lapply(x,'[[', 2) })
    # Combine the braekpoints and distances in one dataframe
    result.dt <- Map(function(bp, dist){
      names(bp) <- names_seq
      names(dist) <- names_seq
      dtlist <- Map(function(per_chr_bp, per_chr_dist){
        data.table(per_chr_bp, per_chr_dist)
      }, bp, dist)
      per_cell_dt <- rbindlist(dtlist, idcol = 'Chr')
      per_cell_dt <- per_cell_dt[order(-per_chr_dist)]
    }, breakpoints, distances)
    
    # Discard irrelevant breakpoints
    pruned_result.dt <- lapply(result.dt, function(x){
      threshold_dist_values(x)
    })
    message("Successfully discarded irrelevant breakpoints")
    
    clusters_pruned <- as.data.table(Map(function(seq_data, bp){
      # Split the readcounts into the different chr
      per_chrom_seq_data <- split(seq_data, peaks$seqnames)
      # Split the calculated breakpoints into the different chr
      per_chrom_bp <- split(bp, bp$Chr)
      # Identify the differnt segments per chromosome after the pruning
      clusters_per_chrom <- Map(function(seq_data2, bp2){
        if(is.null(bp2)){
          return(rep(1, length(seq_data2)))
        } else {
          bp_to_cluster <- sort(c(1,length(seq_data2)+1,bp2$per_chr_bp))
          return(rep(1:length(diff(bp_to_cluster)),diff(bp_to_cluster)))
        }
      }, per_chrom_seq_data[names_seq], per_chrom_bp[names_seq])
      
      # Identify the segments at the genome level after the pruning
      cl <- 0 # counter
      clusters <- vector()
      for(item in clusters_per_chrom){
        item <- item + cl
        clusters <- append(clusters, item)
        cl <- cl + length(unique(item))
      }
      return(clusters)
    },  peaks[, .SD, .SDcols = patterns("cell-")], pruned_result.dt))

    # Assign copy number states to the different "clusters"/segments identified
    somies_ad <- Map(function(seq_data,cluster) {
      assign_gainloss(seq_data, cluster, uq=uq, lq=lq)
    }, peaks[, .SD, .SDcols = patterns("cell-")], clusters_pruned)
    message("Successfully assigned gain-loss")
    saveRDS(somies_ad, file.path(outdir, "cnv_calls.rds"))
  }
  
  somies_ad <- readRDS(file.path(outdir,"cnv_calls.rds"))
  # Write the results to disk
  write_somies.dt <- as.data.table(somies_ad)
  write_somies.dt <- as.data.table(lapply(write_somies.dt, function(x) {
    # x <- factor(x, levels = c(0,1,2), labels=c("Loss", "Normal", "Gain"))
    x <- factor(x, levels = c(0,1,2))
    return(x)
    }))
  write_somies.dt <- as.data.table(cbind(seq=peaks$seqnames, start=peaks$start, 
                                         end=peaks$end, write_somies.dt))
  write.table(write_somies.dt, file = file.path(outdir, "results_table.tsv"), quote = FALSE)
  message("A .tsv file with the results has been written to disk. 
          It contains the copy number states for each cell per bin.
          0 denotes 'Loss', 1 denotes 'Normal', 2 denotes 'Gain'.")
  
  # ----------------------------------------------------------------------------
  # Plotting the result karyogram
  # ----------------------------------------------------------------------------
  
  if(plotKaryo){
    if(is.null(title_karyo)){
      title_karyo <- basename(outdir)
    }
    # Call plotting function to plot and save karyogram
    plot_karyo_gainloss(somies_ad = somies_ad, outdir = outdir, peaks = peaks, uq, lq, title_karyo)
    message("Successfully plotted karyogram")
  }
}
