elbow_finder <- function(dist_values) {
  if(which.max(dist_values)==1){
    y_values <- dist_values
    add_to_x <- 0
  }else {
    split_at_max <- splitAt(dist_values, which.max(dist_values))
    # y_values <- sort(split_at_max[[2]], decreasing = TRUE)
    y_values <- split_at_max[[2]]
    add_to_x <- length(split_at_max[[1]])
  }
  x_values <- c(1:length(y_values))
  max_x_x <- max(x_values)
  max_x_y <- y_values[which.max(x_values)]
  max_y_y <- max(y_values)
  max_y_x <- x_values[which.max(y_values)]
  max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))

  # Creating straight line between the max values
  fit <- stats::lm(max_df$y ~ max_df$x)

  # Distance from point to line
  distances <- c()
  for(i in 1:length(x_values)) {
    distances <- c(distances, abs(stats::coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
  }

  # Max distance point
  x_max_dist <- x_values[which.max(distances)]
  y_max_dist <- y_values[which.max(distances)]

  return(x_max_dist + add_to_x)
}
