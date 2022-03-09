#' @export
#' @export
#' @param somies_ad A list containing the somy per bin of each cell
#' @param outdir Directory where the output karyogram is to be saved
#' @param peaks Dataframe containing bin information
plot_karyo_gainloss <- function(somies_ad, outdir, peaks, uq=NULL, lq=NULL, title_karyo=NULL){
  qc_dt <- data.table()
  qc_dt$spikiness <- sapply(peaks[, .SD, .SDcols = patterns("cell-")], qc.spikiness)
  qc_dt$entropy <- sapply(peaks[, .SD, .SDcols = patterns("cell-")], qc.entropy)
  qc_dt$sumsquares <- unlist(Map(function(counts,somies) {
    qc.sos(counts,somies)
  }, peaks[, .SD, .SDcols = patterns("cell-")], somies_ad))
  qc_dt$libsize <- sapply(peaks[, .SD, .SDcols = patterns("cell-")], sum)
  print(dim(qc_dt))
  somies.dt <- as.data.table(somies_ad)
  write.table(somies.dt, file = file.path(outdir, "results_table.tsv"), quote = FALSE)
  print(dim(somies.dt))
  # somies.dt <- as.data.table(lapply(somies.dt, function(x) {scale(x, center=TRUE, scale=TRUE)}))
  qc_dt$name <- colnames(somies.dt)
  somies.dt$seqnames <- peaks$seqnames
  somies.dt$rn <- as.numeric(rownames(somies.dt))
  somies_melted <- melt(somies.dt, id.vars=c('rn','seqnames'))
  somies_melted$value <- as.factor(paste0(somies_melted$value,'-somy'))
  counts_t <- t(somies.dt[ ,.SD, .SDcols=patterns('cell-')])
  if(nrow(counts_t)>1){
    dist_matrix <- dist(counts_t)
    dist_matrix[is.na(dist_matrix)] <- 0
    hc_counts <- hclust(dist_matrix, method = "ward.D")
    ord <- hc_counts$order
    dhc <- stats::as.dendrogram(hc_counts)
    ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
    ggdndr <- ggplot(ddata$segments) + geom_segment(aes_string(x='x', xend='xend', y='y', yend='yend')) + scale_y_reverse()
    ggdndr <- ggdndr + coord_flip()
    ggdndr <- ggdndr + theme(panel.background=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank(), axis.title=element_blank())
    somies_melted$variable <- factor(somies_melted$variable,
                                     levels = names(somies_ad)[ord])
  }
  somycolours <- c(`0-somy` = "darkorchid3",
                   `1-somy` = "springgreen2",
                   `2-somy` = "red3")
  # text_subtitle <- paste0("Segment lq: ", lq, "Segment uq: ", uq)
  ggsomy <- ggplot(somies_melted, aes(x=rn, y=variable, fill=value)) + geom_tile() +
    facet_grid(cols=vars(seqnames), scales = 'free_x', space = 'free') +
    labs(x="Position in chromosome", fill='Somy', title = title_karyo) +
    scale_fill_manual(values=somycolours) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = 'none',
          strip.text.x = element_text(size = 14),
          axis.title.y = element_blank(),
          axis.text.y = element_blank())

  karyoname <- paste0("Karyogram.png")
  outkaryo <- file.path(outdir, karyoname)

  ggsave(outkaryo, ggsomy, width = 30, height=20, units = "in")
}
