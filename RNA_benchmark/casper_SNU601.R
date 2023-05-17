# ------------------------------------------------------------------------------
# Run CaSpER on the SNU601 cell line
# ------------------------------------------------------------------------------

library(CaSpER)
library(data.table)
library(Seurat)

path_matrix<-"results/inferCNV/input_SNU601_withref/"
path_results<-"results/CaSpER/"

sample_name_run<-"SNU601"

#Load count matrix
count_matrix<-fread(paste0(path_matrix,
                           "count_matrix.txt"))

#Transform it into a numeric matrix
gene_names<-count_matrix$V1
count_matrix$V1<-NULL
count_matrix<-as.matrix(count_matrix)
rownames(count_matrix)<-gene_names

#Get cytoband information
data(hg19_cytoband)

#Get annotation (works with ensembl_gene_id and hgnc_symbol)
annotation <- generateAnnotation(id_type="hgnc_symbol",
                                 genes=gene_names,
                                 ishg19=T, centromere)

#Load metadata
meta_data <- fread(paste0(path_matrix,"sample_annotation.txt"),header=FALSE)
colnames(meta_data)<-c("cell","sample")

#Load count matrix into seurat for filtering and normalization
count_seurat <- CreateSeuratObject(counts = count_matrix, project = "bcc", 
                                   min.cells = 3, min.features = 200)
count_seurat <- NormalizeData(count_seurat , scale.factor = 1e6, 
                              normalization.method = "RC")

log.ge <- as.matrix(count_seurat@assays$RNA@data)
log.ge <- log2(log.ge +1)

rm(count_seurat)
gc()

#Filter gene matrix to contain only annotated cells in the right order)
annotation<-annotation[annotation$Gene %in% rownames(log.ge),]
log.ge <- log.ge[match(annotation$Gene,rownames(log.ge)) , ]
all(rownames(log.ge)==annotation$Gene)

#Extract barcodes of healthy cells as controls
control_cells <- meta_data$cell[meta_data$sample != "SNU601"]

#Load BAFExtract output file
loh <- readBAFExtractOutput(path=paste0(path_results,"SNU601_BAFExtract_af"), 
                            sequencing.type="single-cell",
                            suffix="af")

#Rename files to match sample name
names(loh)<-gsub("_BAFExtract.af","",names(loh))

#Create sample - barcode mapping for LOH annotation
loh_name_mapping <- data.frame(loh.name= meta_data$sample , 
                               sample.name=meta_data$cell)

#Initialize CaSpER object
object <- CreateCasperObject(raw.data=log.ge,
                             loh.name.mapping=loh_name_mapping, 
                             sequencing.type="single-cell", 
                             cnv.scale=3, loh.scale=3, 
                             expr.cutoff=0.1, filter="median", 
                             matrix.type="normalized",
                             annotation=annotation, method="iterative", 
                             loh=loh, 
                             control.sample.ids=control_cells, 
                             cytoband=cytoband)

# Run CaSpER
final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, 
                           method="iterative")

#Save results
saveRDS(final.objects,file=paste0(path_results,"results_",
                                  sample_name_run,"/casper_",
                                  sample_name_run,"_results.rds"))

# summarize large scale events
finalChrMat <- extractLargeScaleEvents(final.objects, thr=0.75)
saveRDS(finalChrMat,file=paste0(path_results,
                                "results_",sample_name_run,"/large_scale_events_",
                                sample_name_run,".rds"))
print("Large scale events extracted")

## plot large scale events
plot.data <- reshape2::melt(finalChrMat)
plot.data$value2 <- "neutral"
plot.data$value2[plot.data$value > 0] <- "amplification"
plot.data$value2[plot.data$value < 0] <- "deletion"
plot.data$value2 <- factor(plot.data$value2, levels = c("amplification",
                                                        "deletion", "neutral"))
plot.data$Var2 <- factor(plot.data$Var2, levels = colnames(finalChrMat))
p <- ggplot(aes(x = Var2, y = Var1, fill = value2), data = plot.data) +
  geom_tile(colour = "white", size = 0.01)+
  labs(x = "",
       y = "") + scale_fill_manual(values = c(amplification = muted("red"),
                                              deletion = muted("blue"),
                                              neutral = "white"))+
  theme_grey(base_size = 6) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.title = element_blank(), strip.text.x = element_blank(),
        legend.text = element_text(colour = "black", size = 7,
                                   face = "bold"),
        legend.key.height = grid::unit(0.8, "cm"), legend.key.width = grid::unit(0.5, "cm"),
        axis.text.x = element_text(size = 5, colour = "black",
                                   angle = -45, hjust = 0),
        axis.text.y = element_text(size = 6, vjust = 0.2, colour = "black"),
        axis.ticks = element_line(size = 0.4),
        plot.title = element_text(colour = "black", hjust = 0,
                                  size = 6, face = "bold"))
ggsave(p,file=paste0(path_results,"large_scale_events_",
                     sample_name_run,".png"),
       width=12,height=10)

#Heatmap with normalized expression
obj <- final.objects[[9]]
plotHeatmap10x(object=obj, fileName=paste0(path_results,
                                           "results_",sample_name_run,
                                           "/heatmap_",sample_name_run,".png"),
             cnv.scale= 3, cluster_cols = F,
             cluster_rows = T, show_rownames = F, only_soi = T)
print("Heatmap saved")

#Segment summary
gamma<-7
segment.summary <- extractSegmentSummary (final.objects)
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
#loh <- segment.summary$all.summary.loh #ignore loss of heterozygosity for the moment

loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.segments<-rbind(loss.final,gain.final)
write.table(all.segments,file=paste0(path_results,
                                     "results_",sample_name_run,
                                     "/segment_summary_",
                                     sample_name_run,".tsv"))
print("Saved segmentation results detailed!")


#Filter the segments for the cancer cells and create a Grange
cancer_cells <- meta_data$cell[meta_data$sample %in% "SNU601"]
segments <- segments[segments$ID %in% cancer_cells,]
seg_grange <-  GRanges(seqnames = Rle(paste0("chr",gsub("p|q", "", segments$seqnames))), 
                       IRanges(segments$start, segments$end)) 

#Create a Grange from the gene annotations
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

saveRDS(ann_gr,file=paste0(path_results,"casper_pseudobulk_",
                           sample_name_run,".rds"))


## plot BAF deviation
plotBAFAllSamples (loh = final.objects[[9]]@loh.median.filtered.data,  
                   fileName=paste0(path_results,
                                   "results_",sample_name_run,
                                   "/LOH_",sample_name_run,".png"))
print("Saved LOH median filtered!")
