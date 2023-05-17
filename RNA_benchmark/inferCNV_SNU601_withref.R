# ------------------------------------------------------------------------------
# Run inferCNV on the SNU601 cell line
# ------------------------------------------------------------------------------

library(infercnv)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="data/input_SNU601/count_matrix.txt",
                                    annotations_file="data/input_SNU601/sample_annotation.txt",
                                    delim="\t",
                                    gene_order_file="data/annotations/gencode_v19_gene_pos.txt",
                                    ref_group_names=c("Patient25","Patient26","Patient27","Patient28","Patient29"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, 
                             out_dir="results/inferCNV/cnv_results_SNU601_withref", 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
			     HMM_type="i6",
			     analysis_mode='subclusters')

