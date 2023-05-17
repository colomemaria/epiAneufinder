# ------------------------------------------------------------------------------
# Run copykat on the SNU601 cell line with reference dataset
# ------------------------------------------------------------------------------

library(copykat)
library(data.table)

data_path<-"data/input_SNU601/"
  
#Load dataset as prepared for inferCNV
data_matrix<-fread(paste0(data_path,"count_matrix.txt"))
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
                        norm.cell.names = annotation$V1[annotation$V2!="SNU601"] , #option to add reference cells
                        sam.name="SNU601_withref", #sample name used for output files
                        distance="euclidean",
                        output.seg="FALSE",
                        plot.genes="TRUE",
                        genome = "hg20",
                        n.cores=1)
