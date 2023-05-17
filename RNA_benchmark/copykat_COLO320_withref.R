# ------------------------------------------------------------------------------
# Run copykat on the COLO320 cell line (all 8 replicates together)
# ------------------------------------------------------------------------------

library(copykat)
library(data.table)

#Define the analyzed replicate
rep_number<-"rep1"
print(paste("Analyzing replicate:",rep_number))

data_path<-"data/colo320_HSR_matrix/"

data_matrix<-fread(paste0(data_path,"colo320hsr_",rep_number,"_colon_ref.tsv.gz"))
#Format into a matrix
gene_names<-data_matrix$V1
data_matrix$V1<-NULL
data_matrix<-as.matrix(data_matrix)
rownames(data_matrix)<-gene_names

#Extract reference cells
annotation<-fread(paste0(data_path,"sample_annotation_colo320hsr_",rep_number,".txt"),
                  header=FALSE)

copykat.test <- copykat(rawmat=data_matrix, #2d matrix with gene expression counts
                        id.type="S", #gene id type (symbol or ensemble)
                        cell.line="no", #if data is from pure cell line
                        ngene.chr=5,
                        LOW.DR = 0.05,
                        UP.DR = 0.1,
                        win.size = 25, #minimal window size for segmentation
                        norm.cell.names = annotation$V1[annotation$V2=="Epithelial"], #specifying reference
                        sam.name=paste0("COLO320_withref_",rep_number), #sample name used for output files
                        distance="euclidean",
                        output.seg="FALSE",
                        plot.genes="TRUE",
                        genome="hg20",
                        n.cores=1)
