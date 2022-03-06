#' @export
threshold_dist_values <- function(result.dt) {
  result.dt$zscores <- scale(result.dt$per_chr_dist, center = TRUE, scale = TRUE)
  result.dt <- result.dt[zscores>0,]
  result.dt$zscores <- NULL
  return(result.dt)
}

#' @export
#' @param seq_data Sequential data - Counts per bin
#' @param cluster Vector showing segment identity
#' @param uq Upper quantile to trim to calculate the cluster means
#' @param lq Lower quantile to trim to calculate the cluster means
assign_gainloss <- function(seq_data, cluster, uq=0.8, lq=0.1) {
  counts.normal <- seq_data / mean(seq_data)
  counts.normal[counts.normal< 0] <- 0
  qus_global <- quantile(seq_data, c(0.01, 0.98))
  cnmean <- sapply(split(counts.normal,cluster), function(x) {
    qus <- quantile(x, c(lq, uq))
    y <- x[x >= qus[1] & x <= qus[2] & x >= qus_global[1] & x <= qus_global[2]]
    if(sum(y) == 0 | length(y) == 0)
      y <- x
    mean(y)
  })
  cnmean_significance <- dplyr::between(scale(cnmean), -1, 1)
  cnmean[cnmean_significance] <- mean(cnmean)
  cnmean.scaled <- cnmean/mean(cnmean)
  cnmean.scaled[cnmean.scaled > 2] <- 2
  if(min(cnmean.scaled) < 0){
    stop()
  }
  CN.states <- round(cnmean.scaled[as.character(cluster)])
  return(CN.states)
}

