# Script to run Copy-scAT on aGBM sample 4349

library(ggplot2)
library(patchwork)
library(CopyscAT)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)

initialiseEnvironment(genomeFile="/work/project/ladcol_005/genomes/hg38_chrom_sizes.tsv",
                      cytobandFile="/work/project/ladcol_005/genomes/hg38_1e+06_cytoband_densities_granges.tsv",
                      cpgFile="/work/project/ladcol_005/genomes/hg38_1e+06_cpg_densities.tsv",
                      binSize=1e6,
                      minFrags=1e4,
                      cellSuffix=c("-1","-2"),
                      lowerTrim=0.5,
                      upperTrim=0.8)

setOutputFile("/work/project/ladcol_005/epiAneufinder/results/CopyscAT/","aGBM_4349_1Mb_5000")

scData <- readInputTable("/work/project/ladcol_005/epiAneufinder/data/ATAC/aGBM/CopyscAT/aGBM_4349_1Mb_5000_copyscat_hg38.tsv")
scData_k_norm <- normalizeMatrixN(scData,logNorm = FALSE,maxZero=2000,imputeZeros = FALSE,blacklistProp = 0.8,blacklistCutoff=125,dividingFactor=1,upperFilterQuantile = 0.95)

#saveRDS(scData_k_norm, "/work/project/ladcol_005/epiAneufinder/data/ATAC/aGBM/aGBM_35_100kb_20000_scData_k_norm_copyscat.rds")

summaryFunction<-cutAverage
scData_collapse<-collapseChrom3N(scData_k_norm,summaryFunction=summaryFunction,binExpand = 1,minimumChromValue = 100,logTrans = FALSE,tssEnrich = 1,logBase=2,minCPG=300,powVal=0.73)

scData_collapse<-filterCells(scData_collapse,minimumSegments = 40,minDensity = 0.1)

graphCNVDistribution(scData_collapse,outputSuffix = "aGBM_violins_1Mb_5000")

median_iqr <- computeCenters(scData_collapse,summaryFunction=summaryFunction)

#identify chromosome-level amplifications
candidate_cnvs<-identifyCNVClusters(scData_collapse,median_iqr,useDummyCells = TRUE,propDummy=0.25,minMix=0.01,deltaMean = 0.03,deltaBIC2 = 0.25,bicMinimum = 0.1, subsetSize=50,fakeCellSD = 0.08, uncertaintyCutoff = 0.55,summaryFunction=summaryFunction,maxClust = 4,mergeCutoff = 3,IQRCutoff= 0.2,medianQuantileCutoff = 0.4)
#cleanup step
candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.5) #= 1.5)
#final results and annotation
final_cnv_list<-annotateCNV4(candidate_cnvs_clean, saveOutput=TRUE,outputSuffix = "clean_cnv",sdCNV = 0.5,filterResults=TRUE,filterRange=0.8)

#smoothedCNVList<-smoothClusters(scDataSampClusters,inputCNVList = final_cnv_list[[3]],percentPositive = 0.4,removeEmpty = FALSE)


library(compiler)
dmRead<-cmpfun(identifyDoubleMinutes)
#
#minThreshold is a time-saving option that doesn't call changepoints on any cell with a maximum Z score less than 4 - you can adjust this to adjust sensitivity of double minute calls (note - lower value = slower)
dm_candidates<-dmRead(scData_k_norm,minCells=10,qualityCutoff2 = 100,minThreshold = 4)

write.table(x=dm_candidates,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"samp_dm.csv"),quote=FALSE,row.names = FALSE,sep=",")

