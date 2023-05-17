# ------------------------------------------------------------------------------
# Run inferCNV on the COLO320 dataset with an external reference (gut cell atlas)
# ------------------------------------------------------------------------------

library(infercnv)

#Specify the directories
data_dir<-"data/colo320_HSR_matrix/"
gene_file<-"data/annotations/hg38_gencode_v27.txt"
out_dir<-"results/output_COLO320/infercnv/"

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(data_dir,"colo320hsr_rep1_colon_ref.tsv.gz"),
                                    annotations_file=paste0(data_dir,"sample_annotation_colo320hsr_rep1.txt"),
                                    delim="\t",
                                    gene_order_file=gene_file,
                                    ref_group_names=c("Epithelial"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir=out_dir, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
                  			     HMM_type="i6",
                  			     analysis_mode='subclusters')

