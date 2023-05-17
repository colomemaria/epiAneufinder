# ------------------------------------------------------------------------------
# Compare copyKat with EpiAneufinder per cell (instead of using the pseudobulk)
# for the COLO320 dataset (replicate 1)
# ------------------------------------------------------------------------------

library(data.table)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)

theme_set(theme_bw())

dataset_name<-"COLO320"
replicate<-"rep1"

# ------------------------------------------------------------------------------
# Datasets - file paths
# ------------------------------------------------------------------------------

#Sample annotation and ref groups of the COLO dataset
sample_file<-"data/colo320_HSR_matrix/input_COLO320/sample_annotation.txt"
ref_groups_file<-"data/colo320_HSR_matrix/input_COLO320/ref_groups.txt"

#Path to EpiAneufinder and WGS results (in same file)
epianeufinder_path<-paste0("results/epiAneufinder_results/COLO320/",
                           "GSM4861367_COLO320HSR_rep1_atac/",
                           "colo320HSP_rep1_results_table.tsv")

#Path to result files for copyKat
copykat_path<-paste0("results/output_COLO320/copykat/",
                     "COLO320_copykat_CNA_raw_results_gene_by_cell.txt")

# ------------------------------------------------------------------------------
# Load datasets
# ------------------------------------------------------------------------------

#Read epiAneufinder results
epianeufinder<-fread(epianeufinder_path)

#Create Grange object (to calculate overlaps downstream)
epianeu_range<-makeGRangesFromDataFrame(epianeufinder,seqnames.field="seq")

#Save cells in the respective matrix
epianeu_cells<-epianeufinder[,5:ncol(epianeufinder)]
epianeu_cells<-as.matrix(epianeu_cells)

#Process copyKat results
copykat_res<-fread(copykat_path)

#Generate a grange object
copykat_res$chromosome_name<-paste0("chr",copykat_res$chromosome_name)
copykat_grange<-makeGRangesFromDataFrame(copykat_res,
                                         seqnames.field="chromosome_name",
                                         start.field="start_position",
                                         end.field="end_position")

#Filter for the cells
copykat_cells<-copykat_res[,8:ncol(copykat_res)]
copykat_cells<-as.matrix(copykat_cells)

#Filter for COLO320 cells
annot<-read.table(sample_file)
colo_cells<-intersect(annot$V1[annot$V2=="COLO320"],colnames(copykat_cells))
copykat_cells<-copykat_cells[,colo_cells]

# ------------------------------------------------------------------------------
# Match copyKat and epiAneufinder cells
# ------------------------------------------------------------------------------

#Match copyKat and epiAneufinder cells
colnames(epianeu_cells)<-gsub("cell-","",colnames(epianeu_cells))

#To check the intersect between both
common_cells<-intersect(colnames(epianeu_cells),colnames(copykat_cells))
print(paste("Cells analyzed by both filtered (survived all QC criteria):",
            length(common_cells)))

#Keep only cells present in both
copykat_cells<-copykat_cells[,common_cells]
epianeu_cells<-epianeu_cells[,common_cells]

#Check overlaps (to get the bin sizes to a similar shape)
overlaps<-as.data.frame(findOverlaps(epianeu_range,copykat_grange,type="any"))

#Take all EpiAneufinder regions overlap with at least one gene
filtered_range<-epianeu_range[unique(overlaps$queryHits),]
epianeu_cells<-epianeu_cells[unique(overlaps$queryHits),]

print(paste("Intersect of genomic regions to evaluate:",
            nrow(epianeu_cells)))

#Average across counts for copykat
copykat_condensed <- sapply(unique(overlaps$queryHits),
                            function(i) {indices<-overlaps$subjectHits[overlaps$queryHits==i]
                            if(length(indices)>1){
                              colMeans(copykat_cells[indices,])
                            } else {
                              copykat_cells[indices,]
                            }})
copykat_condensed<-t(copykat_condensed)

dim(copykat_condensed)

# ------------------------------------------------------------------------------
# Calculate and plot per cell correlation
# ------------------------------------------------------------------------------

#Explore the pearson correlation
cor_per_cell<-sapply(1:ncol(epianeu_cells),function(i)cor(copykat_condensed[,i],
                                                          epianeu_cells[,i]))

g<-ggplot(data.frame(cor=cor_per_cell),aes(x=cor))+geom_histogram(bins=50)+
  xlab("Pearson correlation per cell")+ylab("Number of cells")

ggsave(g, file=paste0("results/plots/",dataset_name,"_",replicate,"_percell_RNA_comparison_corr.png"),
       width=5,height=4)
