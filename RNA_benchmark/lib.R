# ------------------------------------------------------------------------------
# Collection of help functions for reading output files of different methods
# and comparing the results
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(liftOver))

#Function to combine output from two different functions 
#bin size of the first range is taken and method column of the second range is transferred
combine_range_objects<-function(range_1,range_2,method_colname){
  
  if(! method_colname %in% colnames(elementMetadata(range_2))){
    stop(paste("Specified column name",method_colname,"not found in range 2!"))
  }
  #Check overlaps (to get the bin sizes to a similar shape)
  overlaps<-as.data.frame(findOverlaps(range_1,range_2,type="any"))
  
  #Filter first range (keep dimensions)
  range_1<-range_1[unique(overlaps$queryHits)]
  
  #Average score of second method across bins
  elementMetadata(range_1)[method_colname] <- sapply(unique(overlaps$queryHits),
                                                     function(i) mean(elementMetadata(range_2)[overlaps$subjectHits[overlaps$queryHits==i], method_colname]))                                 
  
  #Return range 1 with the new additional column from range2                                                   
  return(range_1)
}

#Read results from inferCNV (6 state HMM model)
read_infercnv_6state_model<-function(gene_path,matrix_path,filter=NULL,
                                    get_combined_column=FALSE){
  res_inferCNV<-fread(matrix_path)
  
  
  #Transform into a numeric matrix
  gene_names<-res_inferCNV$V1
  res_inferCNV$V1<-NULL
  res_inferCNV<-as.matrix(res_inferCNV)
  rownames(res_inferCNV)<-gene_names
  
  #Filter for a certain sample
  if(! is.null(filter)){
    res_inferCNV<-res_inferCNV[,startsWith(colnames(res_inferCNV),filter)]
  }
  
  #Get average results across all cells
  pseudobulk_res_inferCNV<-data.frame(gene_names,
                                      mean_loss=rowMeans(res_inferCNV == 0 | res_inferCNV == 0.5),
                                      mean_base=rowMeans(res_inferCNV==1),
                                      mean_gain=rowMeans(res_inferCNV == 1.5 | res_inferCNV == 2 | res_inferCNV == 3))
  
  #Read gene position information
  gene_pos <- fread(gene_path)
  colnames(gene_pos)<-c("gene","chr","start","end")
  
  #Filter for genes that were not filtered out
  pseudobulk_res_inferCNV<-merge(pseudobulk_res_inferCNV,gene_pos,
                                 by.x="gene_names",by.y="gene",sort=FALSE)
  
  gene_pos_range<-makeGRangesFromDataFrame(pseudobulk_res_inferCNV, keep.extra.columns=TRUE)
  
  if(get_combined_column){
      gene_pos_range$infercnv_mean<-with(gene_pos_range,mean_loss+(mean_base*2)+mean_gain*3)
      gene_pos_range$mean_loss<-NULL
      gene_pos_range$mean_base<-NULL
      gene_pos_range$mean_gain<-NULL
  }
    
  return(gene_pos_range)
}

#Read results from casper (normalized expression)
read_casper<-function(input_file){
  
  casper_complete<-readRDS(input_file)
  casper_score_merged<-with(casper_complete,mean_loss+(mean_base*2)+mean_gain*3)
  casper_complete$casper<-casper_score_merged
  
  return(casper_complete)
}


#Read results from copyKat (normalized expression)
read_copykat<-function(matrix_path,annot_path,ref_path){
  
  suppressWarnings(copykat_res<-fread(matrix_path))
  
  #Generate a grange object with pseudobulk results
  copykat_res$chromosome_name<-paste0("chr",copykat_res$chromosome_name)
  
  #Replace chr23 by chrX (to align with other methods)
  copykat_res$chromosome_name[copykat_res$chromosome_name=="chr23"]<-"chrX"
  
  copykat_grange<-makeGRangesFromDataFrame(copykat_res,
                                           seqnames.field="chromosome_name",
                                           start.field="start_position",
                                           end.field="end_position")
  
  #Get cell matrix
  copykat_cells<-as.matrix(copykat_res[,8:ncol(copykat_res)])
  
  #Remove reference cells
  annot<-fread(annot_path,header=FALSE)
  ref_groups<-fread(ref_path)
  cancer_cells<-annot$V1[!(annot$V2 %in% ref_groups$ref_groups)]
  copykat_cells<-copykat_cells[,colnames(copykat_cells) %in% cancer_cells]
  
  #Calculate the pseudobulk combined over all cells
  copykat_grange$copykat<-rowMeans(copykat_cells)+2
  
  return(list(copykat_grange,ncol(copykat_cells)))
}

# Read result table from Epianeufinder and get mean loss, base and gain values
# seq_col_name: there exist two different versions of EpiAneufinder output file formats, with the first
#               first column called either "seqnames" or "seq" and respective 3 or 4 annotated columns
read_epianeufinder_results<-function(path,seq_col_name="seqnames",get_combined_column=FALSE){
  epianeufinder<-fread(path)
  
  #Create Grange object (to calculate overlaps downstream)
  bin_range<-makeGRangesFromDataFrame(epianeufinder,seqnames.field=seq_col_name)
  
  #New version of fileformat
  if(seq_col_name=="seqnames"){
      start_col<-4
  }
  else if (seq_col_name=="seq"){
      start_col<-5
  } else {
      stop("File format not recognized")
  }
  
  #Get either directly the mean per cell directly or loss/base/gain separately
  if(get_combined_column){
      bin_range$epianeu_mean<-rowMeans(epianeufinder[,start_col:ncol(epianeufinder)])+1
  } else {
      #Getting mean loss, base and gain for each bin over all cells
      bin_range$epiAneufinder_mean_loss<-rowMeans(epianeufinder[,start_col:ncol(epianeufinder)]==0)
      bin_range$epiAneufinder_mean_base<-rowMeans(epianeufinder[,start_col:ncol(epianeufinder)]==1)
      bin_range$epiAneufinder_mean_gain<-rowMeans(epianeufinder[,start_col:ncol(epianeufinder)]==2)
  }
    
  return(bin_range)
}

# Calculate different correlation metrics for two methods 
# (selected by column names "var_name_one" and "var_name_two" 
# in the same data frame "bin_df")
evaluate_correlation_colnames<-function(bin_df,var_name_one,var_name_two,printres=TRUE){
  
  #Calculate correlation values
  pearson_corr<-round(cor(bin_df[[var_name_one]],bin_df[[var_name_two]],
                          method="pearson"),3)
  spearman_corr<-round(cor(bin_df[[var_name_one]],bin_df[[var_name_two]],
                           method="spearman"),3)
  kendall_corr<-round(cor(bin_df[[var_name_one]],bin_df[[var_name_two]],
                          method="kendall"),3)
  
  #Calculate mean square error
  mse <- round(mean((bin_df[[var_name_one]]-bin_df[[var_name_two]])^2),3)
    
  #Print results
  if(printres){
      print(paste("Comparison between",var_name_one,"and",var_name_two))
      print(paste("Pearson correlation:",pearson_corr))
      print(paste("Spearman correlation:",spearman_corr))
      print(paste("Kendell correlation:",kendall_corr))
      print(paste("Mean square error:",mse))
  }
    
  return(c(pearson_corr,spearman_corr,kendall_corr,mse))
  
}

#Convert genomic coordinates from hg19 to hg38 or back (dependent on path_chain)
lift_over_genomes<-function(range_object, path_chain){
  
  #Load the selected chain
  ch = import.chain(path_chain)
  
  #Perform the liftover and return results
  #Remark: some regions are split into multiple smaller regions!!!
  range_object_lifted<-liftOver(range_object,ch)
  range_object_lifted<-unlist(range_object_lifted)

  return(range_object_lifted)
}

#Read WGS results (preprocessed in python to pseudobulk) and transform into an Grange object
read_wgs_results<-function(path,add_chr_symbol=TRUE){
  
  wgs<-fread(path)
  wgs$V1<-NULL
  
  #Convert chromosome names to match EpiAneufinder and inferCNV results
  if(add_chr_symbol){
      wgs$chr<-paste0("chr",wgs$chr)
  }
  
  #Combine loss, base and gain to one number
  wgs$wgs_mean<-with(wgs,loss_wgs+(base_wgs*2)+gain_wgs*3)
  
  #Create Grange results
  wgs_grange<-makeGRangesFromDataFrame(wgs, keep.extra.columns=TRUE)
  return(wgs_grange)
}

