#' @export
splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))

#' Read directly a count matrix saved in mtx format
#' 
#' @param dirpath Directory path that needs to contains three files in that directory:
#'                matrix.mtx(.gz), barcodes.tsv(.gz) and peaks.bed(.gz)
#' @import Matrix
#' @export
readCountMatrix<-function(dir_path){

  #Read count matrix
  if(file.exists(file.path(dir_path,"matrix.mtx"))){
    counts<-readMM(file.path(dir_path,"matrix.mtx"))
  } else if (file.exists(file.path(dir_path,"matrix.mtx.gz"))){
    counts<-readMM(file.path(dir_path,"matrix.mtx.gz"))
  } else {
    stop(paste("Count matrix not found, needs to be named",
               "count.mtx(.gz)!"))
  }
  
  #Add cell barcodes
  if(file.exists(file.path(dir_path,"barcodes.tsv"))){
    barcodes<-fread(file.path(dir_path,"barcodes.tsv"),header=FALSE)
  } else if (file.exists(file.path(dir_path,"barcodes.tsv.gz"))){
    barcodes<-fread(file.path(dir_path,"barcodes.tsv.gz"),header=FALSE)
  } else {
    stop(paste("Barcode file for the count matrix not found, needs to be named",
               "barcodes.tsv(.gz)!"))
  }
  
  colnames(counts)<-barcodes$V1
  
  #Add region information
  if(file.exists(file.path(dir_path,"peaks.bed"))){
    regions<-fread(file.path(dir_path,"peaks.bed"),header=FALSE)
  } else if(file.exists(file.path(dir_path,"peaks.bed.gz"))){
    regions<-fread(file.path(dir_path,"peaks.bed.gz"),header=FALSE)
  } else {
    stop(paste("Peak region file for the count matrix not found, needs to be named",
               "peaks.bed(.gz)!"))    
  }
  
  rownames(counts)<-paste(regions$V1,regions$V2,regions$V3,sep="-")
  
  #Convert into a dgCmatrix
  counts<-as(counts,"dgCMatrix")
  
  return(counts)
}


# "%ni%" <- function(){ Negate("%in%") }

#' @export
qc.spikiness <- function(counts) {
  if (is.null(counts)) {
    return(NA)
  }
  counts <- as.vector(counts)
  sum.counts <- sum(counts)
  spikiness <- sum(abs(diff(counts))) / sum.counts
  return(spikiness)
}

#' @export
qc.entropy <- function(counts) {
  if (is.null(counts)) {
    return(NA)
  }
  counts <- as.vector(counts)
  total.counts <- sum(counts)
  n <- counts/total.counts
  entropy <- -sum( n * log(n) , na.rm=TRUE)
  return(entropy)
}

#' @export
qc.sos <- function(counts, somies) {
  sum(counts - somies) ^ 2
}

#' @export
stateColors <- function(states=c('zero-inflation', paste0(0:10, '-somy'), 'total')) {
  state.colors <- c("zero-inflation"="gray90", "0-somy"="gray90",
                    "1-somy"="darkorchid3", "2-somy"="springgreen2",
                    "3-somy"="red3", "4-somy"="gold2", "5-somy"="navy",
                    "6-somy"="lemonchiffon", "7-somy"="dodgerblue",
                    "8-somy"="chartreuse4", "9-somy"="lightcoral",
                    "10-somy"="aquamarine2", "total"="black")
  states.with.color <- intersect(states, names(state.colors))
  cols <- rep('black', length(states))
  names(cols) <- states
  cols[states.with.color] <- state.colors[states.with.color]
  return(cols)
}
