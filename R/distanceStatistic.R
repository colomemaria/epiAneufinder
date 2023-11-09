#' Function to calculate the AD statistic between two distributions
#' @param x Vector with reads counts left of the breakpoint
#' @param y Vector with reads counts right of the breakpoint
#' @param test Test statistics (either AD, KS or Bhattacharya)
#' @export
dist_ad <- function(x, y, test='AD'){
  if(test=='AD'){
    x <- as.numeric(x); y <- as.numeric(y)
    n <- as.numeric(length(x))
    m <- as.numeric(length(y))
    poolsize <- as.numeric(m+n)
    poolvec <- c(x,y)
    pooldistinct <- unique(sort(poolvec))
    sum_x <- as.numeric(0)
    sum_y <- as.numeric(0)
    for(j in 1:(length(pooldistinct)-1)){
      lj <- as.numeric(sum(x == pooldistinct[j]) + sum(y == pooldistinct[j]))
      mxj <- as.numeric(sum(x <= pooldistinct[j]))
      myj <- as.numeric(sum(y <= pooldistinct[j]))
      bj <- as.numeric(sum(x <= pooldistinct[j]) + sum(y <= pooldistinct[j]))
      denom <- as.numeric(poolsize * (bj * (poolsize - bj)))
      num_x <- lj * (((poolsize * mxj) - (n * bj))^2)
      num_y <- lj * (((poolsize * myj) - (m * bj))^2)
      sum_x <- sum_x + (num_x / denom)
      sum_y <- sum_y + (num_y / denom)
    }
    stat_ad <- (sum_x/n) + (sum_y/m)
    return(stat_ad)
  }
  else if(test=='KS'){
    return(stats::ks.test(x,y)$statistic)
  }
  else if(test=='Bhattacharya'){
    mux <- mean(x)
    muy <- mean(y)
    covx <- cov(as.matrix(x))
    covy <- cov(as.matrix(y))
    return(fpc::bhattacharyya.dist(mux, muy, covx, covy))
  }
  else {
    stop("Please choose one of 'AD', 'KS' or 'Bhattacharya' tests")
  }
}

#' Function to calculate the breakpoints with AD stat given a series of data points
#' 
#' @param seq_data Normalized read counts per window
#' @param minsize Integer. Resolution at the level of ins. Default: 1. Setting it to higher numbers runs the algorithm faster at the cost of resolution
#' @param test Test statistics (either AD, KS or Bhattacharya)
#' @export
seq_dist_ad <- function(seq_data, minsize=3, test='AD') {
  bp1 <- seq(from = 1, to = length(seq_data), by = minsize)
  distlist <- vector()
  for(i1 in 1:length(bp1)){
    data_p <- seq_data[1:bp1[i1]]
    data_q <- seq_data[bp1[i1]:length(seq_data)]
    distlist[i1] <- dist_ad(data_p, data_q, test)
  }
  distlist[is.nan(distlist)] <- 0
  return(distlist)
}
