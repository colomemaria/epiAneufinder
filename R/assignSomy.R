#' Prunes breakpoints
#' 
#' @param results result.dt A data table combining the breakpoints and distances
#' @return A data table after pruning the noisy breakpoints
#' @export
threshold_dist_values <- function(result.dt) {
  result.dt$zscores <- scale(result.dt$per_chr_dist, center = TRUE, scale = TRUE)
  result.dt <- result.dt[zscores>0,]
  result.dt$zscores <- NULL
  return(result.dt)
}

#' Does not assign the full somies but stops before discretizing the states
estimate_foldchange <- function(seq_data, cluster, uq=0.9, lq=0.1, 
                                mean_shrinking=FALSE, trimmed_mean=TRUE) {
  counts.normal <- seq_data / mean(seq_data)
  counts.normal[counts.normal< 0] <- 0
  qus_global <- quantile(seq_data, c(0.01, 0.98))
  # Calculate mean read counts per cluster
  cnmean <- sapply(split(counts.normal,cluster), function(x) {
    qus <- quantile(x, c(lq, uq))
    y <- x[x >= qus[1] & x <= qus[2] & x >= qus_global[1] & x <= qus_global[2]]
    if(sum(y) == 0 | length(y) == 0)
      y <- x
    mean(y)
  })
  
  # Identify clusters/segments with Z scores between -1 and 1
  if(mean_shrinking){
    cnmean_significance <- dplyr::between(scale(cnmean), -1, 1)
    # Collapse these cluster means to the genomic cluster mean to keep multiplicities low
    cnmean[cnmean_significance] <- mean(cnmean)
  }
  
  if(trimmed_mean){
    #cnmean.scaled <- cnmean/mean(cnmean,trim=0.1)
    cnmean.scaled <- cnmean/trimmed_mean_iqr(cnmean)
  } else {
    cnmean.scaled <- cnmean/mean(cnmean)
  }
  

  return(cnmean.scaled[as.character(cluster)])
}

#Filter data based on inter-quantile range (q3-q1) before mean calculation
trimmed_mean_iqr<-function(x){
  q1<-quantile(x,0.25)
  q3<-quantile(x,0.75)
  iqr<-q3-q1
  trimmed_x <- x[x>=(q1-1.5*iqr) & x <=(q3+1.5*iqr)]
  return(mean(trimmed_x))
}

#' Assign CNV state
#' 
#' @param seq_data Sequential data - Counts per bin
#' @param cluster Vector showing segment identity
#' @param uq Upper quantile to trim to calculate the cluster means
#' @param lq Lower quantile to trim to calculate the cluster means
#' @return Copy number states for the different segments
#' @export
assign_gainloss <- function(seq_data, cluster, uq=0.9, lq=0.1, mean_shrinking=FALSE) {
  counts.normal <- seq_data / mean(seq_data)
  counts.normal[counts.normal< 0] <- 0
  qus_global <- quantile(seq_data, c(0.01, 0.98))
  # Calculate mean read counts per cluster
  cnmean <- sapply(split(counts.normal,cluster), function(x) {
    qus <- quantile(x, c(lq, uq))
    y <- x[x >= qus[1] & x <= qus[2] & x >= qus_global[1] & x <= qus_global[2]]
    if(sum(y) == 0 | length(y) == 0)
      y <- x
    mean(y)
  })
  
  # Identify clusters/segments with Z scores between -1 and 1
  if(mean_shrinking){
    cnmean_significance <- dplyr::between(scale(cnmean), -1, 1)
    # Collapse these cluster means to the genomic cluster mean to keep multiplicities low
    cnmean[cnmean_significance] <- mean(cnmean)
  }

  cnmean.scaled <- cnmean/mean(cnmean)
  cnmean.scaled[cnmean.scaled > 2] <- 2
  if(min(cnmean.scaled) < 0){
    stop()
  }
  CN.states <- round(cnmean.scaled[as.character(cluster)])
  return(CN.states)
}



