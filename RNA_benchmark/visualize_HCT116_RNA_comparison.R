# ------------------------------------------------------------------------------
# Compare RNA methods (inferCNV, copyKat) with EpiAneufinder and scWGS
# results for the HCT116 dataset
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)

source("lib.R")

dataset_name<-"HCT116"

# ------------------------------------------------------------------------------
# Datasets - file paths
# ------------------------------------------------------------------------------

#Path to results files for inferCNV
infercnv_path<-"results/inferCNV/"
input_gene_pos_inferCNV <- "data/annotations/gencode_v19_gene_pos.txt"
input_inferCNV<-paste0(infercnv_path,"cnv_results_HCT116_withref/",
                       "infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt")

#Path to results files for CaSpER
casper_path<-paste0("results/casper/results_HCT116/",
                    "HCT116_casper_pseudobulk_aggregate.RDS")

#Path to result files for copyKat
copykat_path<-paste0("results/copyKat/HCT116_withref/",
                     "HCT116_withref_copykat_CNA_raw_results_gene_by_cell.txt")

#Path to EpiAneufinder results
epianeufinder_path<-"results/epiAneufinder_results/HCT116/results_table.tsv"

#Path to WGS results
wgs_path<-"results/HCT_scWGS_results/wgs_results_aneufinder_pseudobulk.csv"

#Liftover necessary for inferCNV
path_chain<-"data/annotations/hg19ToHg38.over.chain"

# ------------------------------------------------------------------------------
# Load datasets
# ------------------------------------------------------------------------------

#Read epianeufinder results
epianeu_grange<-read_epianeufinder_results(epianeufinder_path, seq_col_name="seq",
                                           get_combined_column=TRUE)

#Read scWGS data
wgs_range<-read_wgs_results(wgs_path,add_chr_symbol=FALSE)

combined_range<-combine_range_objects(epianeu_grange,wgs_range,
                                      method_colname="wgs_mean")

#Load grange object with casper pseudobulk results (precalculated because of long runtime)
casper_grange <- read_casper(casper_path)

#Combine EpiAneufinder and CaSpER results
combined_range<-combine_range_objects(combined_range,casper_grange,
                                      method_colname="casper")

#Add inferCNV results
infercnv_grange <- read_infercnv_6state_model(input_gene_pos_inferCNV,input_inferCNV)

#Liftover to from hg19 to hg38
infercnv_grange <- lift_over_genomes(infercnv_grange,path_chain)

#Combine inferCNV results to pseudobulk
infercnv_grange$infercnv<-infercnv_grange$mean_loss + infercnv_grange$mean_base*2 +
  infercnv_grange$mean_gain*3

combined_range<-combine_range_objects(combined_range,infercnv_grange,
                                      method_colname="infercnv")

#Add copyKat results
copykat_res<-fread(copykat_path)

#Generate a grange object with pseudobulk results
copykat_res$chromosome_name<-paste0("chr",copykat_res$chromosome_name)
copykat_grange<-makeGRangesFromDataFrame(copykat_res,
                                         seqnames.field="chromosome_name",
                                         start.field="start_position",
                                         end.field="end_position")

#Filter for HCT cells (remove all reference cells)
copykat_hct_cells<-startsWith(colnames(copykat_res),"scCAT_HCT116_")
copykat_grange$copykat_mean<-rowMeans(copykat_res[,copykat_hct_cells,with=FALSE])

combined_range<-combine_range_objects(combined_range,copykat_grange,
                                      method_colname="copykat_mean")

#Shift copykat results as they intially range from -1 to 1
combined_range$copykat_mean<-combined_range$copykat_mean+2

#Save file with combined datasets
write.table(combined_range,file="results_RNA_comparison/results_RNA_comparisons_HCT116.tsv",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

# ------------------------------------------------------------------------------
# Combined line plot
# ------------------------------------------------------------------------------

#Vector for renaming methods (official published names)
method_names<-setNames(c("EpiAneufinder","scWGS","CaSpER","InferCNV","copyKat"),
                       c("epianeu","wgs","casper","infercnv","copykat"))

#Get results
combined_methods<-elementMetadata(combined_range)

#Rename methods (remove "_mean" and "_results")
colnames(combined_methods)<-gsub("_.*","",colnames(combined_methods))

#Add position information
combined_methods$chr<-factor(combined_range@seqnames,
                             levels=combined_range@seqnames@values)
combined_methods$start_position<-combined_range@ranges@start
combined_methods<-as.data.frame(combined_methods)

#Add an artifical count through the whole genome and 
#get start positions for each new chromosome
combined_methods$counted_pos<-1:nrow(combined_methods)
chr_boundries<-combined_methods%>%
  group_by(chr)%>%
  summarize(start_chr=min(counted_pos),
            mean_chr=mean(counted_pos))

#Scale every dataset to have diploid values at 0 and a standard deviation of 1
scaled_methods<-combined_methods
for(method in names(method_names)){
  scaled_methods[,method]<-(scaled_methods[,method]-2) / 
    sd(scaled_methods[,method])
}

plot_data<-reshape2::melt(scaled_methods,
                          id.vars=c("chr","start_position","counted_pos"))

#Rename variable names
plot_data$variable<-method_names[plot_data$variable]

#Order them respectively
plot_data$variable<-factor(plot_data$variable,
                           levels=c("scWGS","EpiAneufinder","Copy-scAT",
                                    "InferCNV","copyKat","CaSpER"))

method_colors<-c("#2ca02c","#ff7f0e","#1f77b4","#9467bd","#e377c2")

#Combine heatmap and line plot
g.1<-ggplot(plot_data,aes(x=counted_pos,y=value,color=variable))+geom_line()+
  geom_vline(xintercept = chr_boundries$start_chr)+
  ylab("Normalized score")+
  scale_color_manual("Method",values=method_colors)+
  scale_x_continuous(breaks=chr_boundries$mean_chr,
                     labels=chr_boundries$chr)+
  coord_cartesian(xlim=c(1,max(plot_data$counted_pos)),expand=FALSE)+
  theme_bw()+
  theme(legend.position="top",
        #legend.key.size = unit(2,"line"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        text=element_text(size=21))

g.2<-ggplot(plot_data,aes(x=counted_pos,y=variable,fill=value))+geom_tile()+
  theme_bw()+
  scale_fill_gradient2("Score",low = "darkblue",
                        mid = "white",high = "darkred",midpoint = 0)+
  xlab("Chromosome position")+ylab("Method")+
  geom_vline(xintercept = chr_boundries$start_chr)+
  scale_x_continuous(breaks=chr_boundries$mean_chr,
                     labels=chr_boundries$chr)+
  scale_y_discrete(limits=rev)+
  coord_cartesian(xlim=c(1,max(plot_data$counted_pos)),expand=FALSE)+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
        text=element_text(size=21))
g<-ggarrange(g.1,g.2,ncol=1,align="v",heights=c(0.3,0.7))
ggsave(g, file=paste0("results/plots/",dataset_name,"_RNA_comparison.pdf"),
       width=20,height=10)
ggsave(g, file=paste0("results/plots/",dataset_name,"_RNA_comparison.png"),
       width=20,height=10)

# ------------------------------------------------------------------------------
# Correlation calculation and plot
# ------------------------------------------------------------------------------

#Create a heatmap which methods are how similar
methods_ordered<-setNames(c("scWGS","EpiAneufinder","InferCNV","CaSpER","copyKat"),
                          c("wgs","epianeu","infercnv","casper","copykat"))
compare_methods<-NULL
for(m1 in 1:(length(methods_ordered)-1)){
  for(m2 in (m1+1):length(methods_ordered)){
    res<-evaluate_correlation_colnames(combined_methods,
                                       names(methods_ordered)[m1],
                                       names(methods_ordered)[m2],
                                       printres=FALSE)
    compare_methods<-rbind(compare_methods,
                           data.frame(method1=methods_ordered[m1],
                                      method2=methods_ordered[m2],
                                      pearson=res[1],
                                      spearman=res[2]))
  }
  
  #Add the diagonal for easier visualization
  compare_methods<-rbind(compare_methods,
                         data.frame(method1=methods_ordered[m1],
                                    method2=methods_ordered[m1],
                                    pearson=1,
                                    spearman=1))
}

#Add the diagonal for the last element
compare_methods<-rbind(compare_methods,
                       data.frame(method1=methods_ordered[m2],
                                  method2=methods_ordered[m2],
                                  pearson=1,
                                  spearman=1))

compare_methods$method1<-factor(compare_methods$method1,
                                levels=methods_ordered)
compare_methods$method2<-factor(compare_methods$method2,
                                levels=methods_ordered)

g<-ggplot(compare_methods,aes(x=method2,y=method1,fill=pearson))+
  geom_tile()+
  theme_bw()+
  ggtitle(paste(dataset_name,"dataset"))+
  geom_text(aes(label=round(pearson,3),
                color=ifelse(pearson<0.6,'white','black')),size=3)+
  scale_color_manual(values=c("black","white"))+
  scale_y_discrete(limits=rev)+
  xlab("Method")+
  ylab("Method")+
  scale_fill_viridis("Pearson\ncorrelation",limits=c(0,1))+
  guides(color="none")+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))

print(g)
ggsave(g, file=paste0("results/plots/",dataset_name,"_RNA_comparison_corr.pdf"),
       width=6,height=4.5)
