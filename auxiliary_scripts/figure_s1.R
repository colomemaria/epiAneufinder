library(data.table)
library(atacCNV)
library(plyr)
library(ggplot2)
library(cowplot)
library(gridExtra)

outdir <- "HCT116_latest"

counts <- readRDS(file.path(outdir,"count_summary.rds"))
peaks <- as.data.table(assays(counts)$counts)
colnames(peaks) <- paste0('cell-', colnames(peaks))
rowinfo <- as.data.table(rowRanges(counts))
peaks <- cbind(rowinfo, peaks)

corrected_counts <- readRDS(file.path(outdir,"counts_gc_corrected.rds"))
peaks <- cbind(rowinfo, corrected_counts)

zeroes_per_bin <- peaks[, rowSums(.SD==0), .SDcols = patterns("cell-")]
ncells <- length(grep("cell-", colnames(peaks)))
threshold_blacklist_bins <- 0.85
peaks <- peaks[zeroes_per_bin<(threshold_blacklist_bins*ncells)]


somies_ad <- readRDS(file.path(outdir,"cnv_calls.rds"))
somies.dt <- as.data.table(somies_ad)

set.seed(123)
randomcells<- sample(x=names(somies.dt), size=50)
singlecellplot <- list()
for(i1 in randomcells){
  plot.dt <- peaks[, .(chr=seqnames, counts=get(i1))]
  plot.dt$somy <- paste0(somies_ad[[i1]], '-somy')
  plot.dt <- plot.dt[which(plot.dt$counts < quantile(plot.dt$counts, 0.999) &
                             plot.dt$counts > quantile(plot.dt$counts, 0.01))]
  plot.dt$rn <- as.numeric(rownames(plot.dt))
  plot.dt$chr <- gsub("chr","", plot.dt$chr)
  plot.dt$chr <- factor(plot.dt$chr, levels = c(1:22))
  titletext <- paste0(i1, " : Library size - ", sum(plot.dt$counts))
  somycolours <- c(`0-somy` = "darkorchid3",
                   `1-somy` = "springgreen2",
                   `2-somy` = "red3")
  numwindows <- as.data.frame(table(plot.dt$somy))
  colnames(numwindows) <- c("State", "Number of bins")
  numwindows$State <- gsub("0-somy","loss", numwindows$State)
  numwindows$State <- gsub("1-somy","normal", numwindows$State)
  numwindows$State <- gsub("2-somy","gain", numwindows$State)
  
  libsize <- as.data.frame(sum(plot.dt$counts))
  colnames(libsize) <- c("Total fragments")
  
  profileplot <- ggplot(plot.dt, aes(x=rn, y=counts)) + geom_point(alpha=0.2) +
    geom_line(aes(y=zoo::rollapply(counts, 200, mean, fill=NA)), color='red', size=1) +
    facet_grid(cols=vars(chr), scales = "free_x", space = "free_x") +
    labs(y="Counts") +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none')
  somyplot <- ggplot(plot.dt, aes(x=rn, y=counts)) + geom_point(alpha=0.2, aes(color=somy)) +
    facet_grid(cols=vars(chr), scales = "free_x", space = "free_x") +
    labs(y="Counts") +
    scale_color_manual(labels = c("loss", "normal", "gain"), values=somycolours) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none')
  somydensity <- ggplot(plot.dt, aes(counts, fill=somy, colour=somy)) +
    geom_density(alpha=0.7) +
    scale_color_manual(labels = c("loss","normal","gain"), values=somycolours) +
    scale_fill_manual(labels = c("loss","normal","gain"), values=somycolours) +
    labs(x="Counts per bin", y="Density", title = titletext) +
    theme(panel.grid.major = element_blank(),
          panel.background = element_blank(),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = 14))
  
  libsizeplot <- ggplot() + annotation_custom(tableGrob(libsize, rows=NULL)) +
    theme(panel.grid.major = element_blank(),
          panel.background = element_blank(),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = 14))
  
  tableplot <- ggplot() + annotation_custom(tableGrob(numwindows, rows=NULL)) +
    theme(panel.grid.major = element_blank(),
          panel.background = element_blank(),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = 14))
  
  singlecellplot[[i1]] <- cowplot::plot_grid(plotlist = list(cowplot::plot_grid(plotlist = list(somydensity, libsizeplot, tableplot), nrow = 1),
                                                             profileplot, 
                                                             somyplot), 
                                             ncol = 1)
}

pdf("hct.pdf", width = 25, height = 15)
for(scp in singlecellplot){
  grid::grid.draw(scp)
}
dev.off()

