# ------------------------------------------------------------------------------
# Run copykat on the HCT116 cell line (from scCAT dataset) with reference dataset
# ------------------------------------------------------------------------------

library(copykat)
library(data.table)

#Load dataset as prepared for inferCNV
data_path<-"data/HCT_processed_matrix/"
data_matrix<-fread(paste0(data_path,"HCT116_withref.tsv"))
#Format into a matrix
gene_names<-data_matrix$V1
data_matrix$V1<-NULL
data_matrix<-as.matrix(data_matrix)
rownames(data_matrix)<-gene_names

#Extract reference cells
annotation<-fread(paste0(data_path,"sample_annotation.txt"),
                  header=FALSE)

copykat.test <- copykat(rawmat=data_matrix, #2d matrix with gene expression counts
                        id.type="S", #gene id type (symbol or ensemble)
                        cell.line="no", #if data is from pure cell line
                        ngene.chr=5,
                        LOW.DR = 0.05,
                        UP.DR = 0.1,
                        win.size = 25, #minimal window size for segmentation
                        norm.cell.names = annotation$V1[annotation$V2!="HCT116"], 
                        sam.name="HCT116_withref", #sample name used for output files
                        distance="euclidean",
                        output.seg="FALSE",
                        plot.genes="TRUE",
                        genome="hg20",
                        n.cores=1)
