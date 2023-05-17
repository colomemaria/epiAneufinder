# ------------------------------------------------------------------------------
# Run inferCNV on the HCT116 cell line
# ------------------------------------------------------------------------------

library(infercnv)

#Specify the directories
data_dir<-"data/HCT_processed_matrix/"
gene_file<-"data/annotations/gencode_v19_gene_pos.txt"
out_dir<-"results/inferCNV/cnv_results_HCT116_withref"

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(data_dir,"HCT116_withref.tsv"),
                                    annotations_file=paste0(data_dir,"sample_annotation.txt"),
                                    delim="\t",
                                    gene_order_file=gene_file,
                                    ref_group_names=c("ILC","B cell","Myeloid cell","CD8 T cell",
                                                      "CD4 T cell","Epithelial cell","Fibroblast"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir=out_dir, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
                  			     HMM_type="i6",
                  			     analysis_mode='subclusters')

