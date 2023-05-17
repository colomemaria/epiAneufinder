# ------------------------------------------------------------------------------
# Run CaSpER on the HCT116 dataset
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(CaSpER)
library(data.table)
library(Seurat)

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

input_file<-"data/HCT_processed_matrix/HCT116_withref.tsv"
input_annotations<-"data/HCT_processed_matrix/sample_annotation.txt"
input_gene_annot<-"annotation_genes_HCT116.tsv"
input_af<-"results/CaSpER/HCT116_BAFExtract_af/HCT116.af"

input_segment_gamma<-7

output_casper_object<-"results_HCT116/HCT116_casper_object_final.RDS"
output_casper_segments<-"results_HCT116/HCT116_casper_segments.tsv"
output_casper_pseudobulk<-"results_HCT116/HCT116_casper_pseudobulk_aggregate.RDS"
 
output_plot_density<-"results_HCT116/HCT116_casper_densities.pdf"
output_plot_large_events<-"results_HCT116/HCT116_casper_large_events.png"
output_plot_heatmap<-"results_HCT116/HCT116_casper_heatmap.png"
output_plot_baf<-"results_HCT116/HCT116_casper_baf.png"

# ------------------------------------------------------------------------------
print("Execute Casper")
# ------------------------------------------------------------------------------

#Load count matrix
count_matrix<-fread(input_file)

#Transform it into a numeric matrix
gene_names<-count_matrix$V1
count_matrix$V1<-NULL
count_matrix<-as.matrix(count_matrix)
rownames(count_matrix)<-gene_names

#Get cytoband information
data(hg38_cytoband)

#Load metadata
meta_data <- fread(input_annotations,header=FALSE)
colnames(meta_data)<-c("cell","sample")

#Load count matrix into seurat for filtering and normalization
#(following the respective single cell tutorial)
count_seurat <- CreateSeuratObject(counts = count_matrix, project = "bcc",
                                   min.cells = 3, min.features = 200)

count_seurat <- NormalizeData(count_seurat , scale.factor = 1e6,
                              normalization.method = "RC")

log.ge <- as.matrix(count_seurat@assays$RNA@data)
log.ge <- log2(log.ge +1)

rm(count_seurat)
gc()

#Filter gene matrix to contain only annotated cells in the right order
annotation<-read.table(input_gene_annot,
                       header=TRUE,stringsAsFactors = FALSE)
annotation<-annotation[annotation$Gene %in% rownames(log.ge),]
log.ge <- log.ge[match(annotation$Gene,rownames(log.ge)) , ]
all(rownames(log.ge)==annotation$Gene)

# ------------------------------------------------------------------------------
print("Extract CNV segments with score")
# ------------------------------------------------------------------------------

#Load final object again (issue in first run)
final.objects <- readRDS(output_casper_object)

print(paste("Chosen thresold to define segments (larger or equal):",
            input_segment_gamma))
segment.summary <- extractSegmentSummary(final.objects)
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
#loh <- segment.summary$all.summary.loh #ignore loss of heterozygosity for the moment

loss.final <- loss[loss$count>=input_segment_gamma, ]
gain.final <- gain[gain$count>=input_segment_gamma, ]
segments<-rbind(loss.final,gain.final)
write.table(segments,quote=FALSE, sep="\t",file=output_casper_segments)
print("Saved segmentation results detailed!")

# ------------------------------------------------------------------------------
print("Create pseudobulk aggregate over all cells")
# ------------------------------------------------------------------------------

#Filter the segments for the cancer cells and create a Grange
cancer_cells <- meta_data$cell[meta_data$sample %in% "HCT116"]
segments <- segments[segments$ID %in% cancer_cells,]
seg_grange <-  GRanges(seqnames = Rle(paste0("chr",gsub("p|q", "", segments$seqnames))), 
                       IRanges(segments$start, segments$end)) 

#Create a Grange from the gene annotations
annotation$Chr<-paste0("chr",annotation$Chr)
ann_gr <- makeGRangesFromDataFrame(annotation, 
                                   keep.extra.columns = TRUE, seqnames.field="Chr")

#Get gene-wise CNV annotations following the code from the CaSpER tutorial
genes <- splitByOverlap(ann_gr, seg_grange, "GeneSymbol")
genes_ann <- lapply(genes, function(x) x[!(x=="")])
rna_matrix <- gene.matrix(seg=segments, all.genes=unique(annotation$GeneSymbol), 
                          all.samples=cancer_cells, 
                          genes.ann=genes_ann)

#Get average results across all cells
ann_gr$mean_loss<-rowMeans(rna_matrix == -1)
ann_gr$mean_base<-rowMeans(rna_matrix == 0)
ann_gr$mean_gain<-rowMeans(rna_matrix == 1)

saveRDS(ann_gr,file=output_casper_pseudobulk)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
