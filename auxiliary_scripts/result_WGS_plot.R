library(plyranges)
library(data.table)
library(atacCNV)

# Split the bed file into multiple files
# First did a grep for "track" and saw that there should be 15 files 
# csplit -k HCT116_GFP_72h_2N_binsize_1e+06_stepsize_5e+05_refined_curated_CNV.bed /track/ {15}

# These bed files can't be converted to GRanges
# Will use the split bed files
# bedfile_cnv <- "Browser_files_curated_samples_Standard/Browser_files_Edivisive_ccc_95_rbc_5/HCT116_GFP_72h_2N_binsize_1e+06_stepsize_5e+05_refined_curated_CNV.bed"
# bedfile_cnv <- fread(bedfile_cnv)

bedfiles <- list.files(path="Browser_files_curated_samples_Standard/Browser_files_Edivisive_ccc_95_rbc_5/",
           pattern="xx", full.names = TRUE)

# Reading in the windows from atacCNV
outdir <- "snu601atac"
counts <- readRDS(file.path(outdir,"count_summary.rds"))
windows <- rowRanges(counts)
mcols(windows) <- NULL

list_cells = list()
for(bedfile in bedfiles) {
  print(basename(bedfile))
  gr <- read_bed(bedfile)
  gr$score <- NULL
  gr$itemRgb <- NULL
  gr$thick <- NULL
  windowsWGS <- mergeByOverlaps(windows, gr)
  windowsWGS$gr <- NULL
  windowsWGS <- windowsWGS[!duplicated(windowsWGS$windows), ]
  windowsWGS <- as.data.table(windowsWGS)
  windowsWGS$windows.width <- NULL
  windowsWGS$windows.strand <- NULL
  setnames(windowsWGS, "name", basename(bedfile))
  list_cells[[basename(bedfile)]] <- windowsWGS
  print(dim(list_cells[[basename(bedfile)]]))
}

somy_WGS <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c('windows.seqnames', 'windows.start', 'windows.end'), all.x = TRUE),
       list_cells)
somies_WGS.dt <- as.data.table(somy_WGS)
somies_WGS.dt$rn <- as.numeric(rownames(somies_WGS.dt))
somies_WGS_melted <- melt(somies_WGS.dt, id.vars=c('rn','windows.seqnames','windows.start','windows.end'))
somies_WGS_melted$value <- as.factor(somies_WGS_melted$value)
somies_WGS_melted$variable <- factor(somies_WGS_melted$variable,
                                     levels = grep('xx', names(somies_WGS.dt), value=TRUE))
# title_karyo <- "WGS"
# ggsomy_somies <- ggplot(somies_WGS_melted, aes(x=rn, y=variable, fill=value)) + geom_tile() +
#   facet_grid(cols=vars(windows.seqnames), scales = 'free_x', space = 'free') +
#   labs(x="Position in chromosome", fill='Somy', title = title_karyo) +
#   scale_fill_manual(values=stateColors(states = unique(somies_WGS_melted$value))) +
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         legend.position = 'none',
#         strip.text.x = element_text(size = 14),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank())


levels(somies_WGS_melted$value)[levels(somies_WGS_melted$value) == "1-somy"] <- "0-somy"
levels(somies_WGS_melted$value)[levels(somies_WGS_melted$value) == "2-somy"] <- "1-somy"
levels(somies_WGS_melted$value)[levels(somies_WGS_melted$value) != "0-somy" & levels(somies_WGS_melted$value) != "1-somy" ] <- "2-somy"

somycolours <- c(`0-somy` = "darkorchid3",
                 `1-somy` = "springgreen2",
                 `2-somy` = "red3")
title_karyo <- "WGS of HCT116"

ggsomy_states <- ggplot(somies_WGS_melted, aes(x=rn, y=variable, fill=value)) + geom_tile() +
  facet_grid(cols=vars(windows.seqnames), scales = 'free_x', space = 'free') +
  labs(x="Position in chromosome", fill='Somy', title = title_karyo) +
  scale_fill_manual(values=somycolours) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'none',
        strip.text.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
