#' Identify breakpoints
#'
#' @param seq_data Sequential data - Counts per bin for a chromosome
#' @param minsize Integer. Resolution at the level of ins. Default: 1. Setting it to higher numbers runs the algorithm faster at the cost of resolution
#' @param k Integer. Find 2^k segments per chromosome
#' @param test Test statistics
#' @param minsizeCNV Integer. Number of consecutive bins to constitute a possible CNV
#' @export
getbp <- function(seq_data, minsize=1, k=3, test='AD', minsizeCNV=5){
  bp <- vector()
  distbp <- vector()
  for(iter in 1:k) {
    
    #Split the vector on the respective breakpoints into segments
    seq_k_data <- splitAt(seq_data, (bp))
    
    #Calculate the distance for each possible breakpoint within the segment
    dist_vect <- lapply(seq_k_data, function(x) {seq_dist_ad(x, minsize, test)})
    
    if(length(bp)==0){
      add_to_bp <- 0
    } else{
      add_to_bp <- c(0, cumsum(sapply(seq_k_data[-length(seq_k_data)], length)))
    }

    #Identify the positions with the maximal distance within each segment
    bp_per_seq <- sapply(dist_vect, function(x) { which.max(x)*minsize })
  
    bp_per_seq[which(bp_per_seq > sapply(seq_k_data, length))] <- 
      sapply(seq_k_data, length)[which(bp_per_seq > sapply(seq_k_data, length))]
    
    bp_neighbors <-  as.vector(sapply(bp, function(x){seq(x-minsizeCNV,x+minsizeCNV)}))
    # bp_per_seq[((bp_per_seq+add_to_bp) %in% bp_neighbors)] <- NA
    sig_per_seq <- !((bp_per_seq+add_to_bp) %in% bp_neighbors)
    # sig_per_seq <- sapply(dist_vect, max)>10
    
    bp_per_seq <- bp_per_seq + add_to_bp
    dist_per_seq <- sapply(dist_vect, max)
    if(sum(sig_per_seq!=0)){
      bp <- append(bp, bp_per_seq[sig_per_seq])
      # distbp <- append(distbp, sapply(dist_vect[sig_per_seq], max))
      distbp <- append(distbp, dist_per_seq[sig_per_seq])
    } else {
      break()
    }
  }
  
  if (minsizeCNV!=0){
      bp_ends <- c(1:minsizeCNV, seq(length(seq_data)-minsizeCNV, length(seq_data)))
     } else{
      bp_ends = c(1, length(seq_data))
     }
  distbp <- distbp[!(bp %in% bp_ends)]
  bp <- bp[!(bp %in% bp_ends)]
  res.df <- data.frame(bp, distbp)
  res.df <- res.df[order(-res.df$distbp),]
  res <- list()
  res[['bp']] <- res.df$bp
  res[['ad.distance']] <- res.df$distbp
  # res[['estimates']] <- unique(res$bp[1:(elbow_finder(res$ad.distance)-1)])
  bp_to_cluster <- sort(c(1,length(seq_data)+1,res$bp))
  res[['cluster']] <- rep(1:length(diff(bp_to_cluster)),diff(bp_to_cluster))
  # res[['somy']] <- assign_somy(seq_data,res$cluster)
  return(res)
}
