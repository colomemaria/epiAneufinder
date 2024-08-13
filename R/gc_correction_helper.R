#' Help function for GC correction per cell, in case that instead of the loess function, 
#' a standard implemetation is used with a quadratic fit for each polynome
#'
#' @param counts Vector with counts of one cell across all bins
#' @param gc_content Vector with gc content per bin
#' @param intervals_per_bin Vector assigning bin to GC interval
#' @param mean_gc_interval Mean GC content of each interval
#' @export
quadratic_gc_correction<-function(counts,gc_content,intervals_per_bin,
                                  mean_gc_interval){
  
  intervals <- sort(unique(intervals_per_bin))
  mean.counts.global <- mean(counts, trim=0.05)
  correction.factors <- NULL
  weights <- NULL
  
  for (interval in intervals) {
    mask <- intervals_per_bin==interval
    counts.with.same.GC <- counts[mask]
    weights[as.character(mean_gc_interval[interval])] <- length(counts.with.same.GC)
    mean.counts.with.same.GC <- mean(counts.with.same.GC, na.rm=TRUE, trim=0.05)
    if (mean.counts.with.same.GC == 0) {
      correction.factor <- 0
    } else {
      correction.factor <-  mean.counts.global / mean.counts.with.same.GC
    }
    correction.factors[as.character(mean_gc_interval[interval])] <- correction.factor
  }
  
  y <- correction.factors
  x <- as.numeric(names(y))
  w <- weights
  df <- data.frame(x,y,weight=w)
  weight <- w	# dummy assignment to pass R CMD check, doesn't affect the fit
  fit <- stats::lm(y ~ poly(x, 2, raw=TRUE), data=df, weights=weight)
  fitted.correction.factors <- stats::predict(fit, data.frame(x=mean_gc_interval[intervals]))
  names(fitted.correction.factors) <- intervals
  
  corrected_counts<-counts
  for (interval in intervals) {
    mask <- which(intervals_per_bin==interval)
    correction.factor <- fitted.correction.factors[as.character(interval)]
    corrected_counts[mask] <- counts[mask] * correction.factor
  }
  
  return(as.integer(round(corrected_counts)))
  
}