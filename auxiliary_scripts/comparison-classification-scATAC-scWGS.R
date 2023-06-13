
rm(list=ls())
library(plotly)

path="/folder/"
seg <- read.table(paste0(path,"SNU_pseudoprofiles.tsv"), sep = '\t', header = TRUE)
colnames(seg)<-c("chr", "start", "end", "scATAC", "scDNA")

lim_low <- seq(from = 1.0, to = 1.995, by = 0.005) #the values of the CNV is between 1 and 3, between 1 and 2 is a loss
lim_high <- seq(from = 2.0, to = 2.9, by = 0.005) #the values of the CNV is between 1 and 3, between 2 and 3 is a gain
precision_all <- recall_all <- rand <- NULL
for(i in 1:length(lim_low)){
  print(i)
  for(j in 1:length(lim_high)){
    #discretize the data here: above the limit it is a CNV of "gain" (labelled as -3), below a CNV of "loss" (labelled as -1), in between a CNV of "normal" (labelled as -2)
    atac2 <- seg$scATAC
    atac2[which(atac2 < lim_low[i])] <- -1
    atac2[which(atac2 > lim_high[j])] <- -3
    atac2[which(atac2 > 0.0)] <- -2
    wgs2 <- seg$scDNA
    wgs2[which(wgs2 < lim_low[i])] <- -1
    wgs2[which(wgs2 > lim_high[j])] <- -3
    wgs2[which(wgs2 > 0.0)] <- -2
        
    #make a confusion table:      
    index1atac=which(atac2==-1)
    index2atac=which(atac2==-2)
    index3atac=which(atac2==-3)
    #confusionMatrix
    cM <- matrix(nrow=3, ncol=3)
    cM[1,1] = length(which(wgs2[index1atac]==-1))
    cM[1,2] = length(which(wgs2[index1atac]==-2))
    cM[1,3] = length(which(wgs2[index1atac]==-3))
        
    cM[2,1] = length(which(wgs2[index2atac]==-1))
    cM[2,2] = length(which(wgs2[index2atac]==-2))
    cM[2,3]= length(which(wgs2[index2atac]==-3))
        
    cM[3,1] = length(which(wgs2[index3atac]==-1))
    cM[3,2] = length(which(wgs2[index3atac]==-2))
    cM[3,3] = length(which(wgs2[index3atac]==-3))
        
    #precision
    pre = c(cM[1,1]/rowSums(cM)[1], cM[2,2]/rowSums(cM)[2], cM[3,3]/rowSums(cM)[3])
    precision_all <- rbind(precision_all, c(lim_low[i], lim_high[j],lim_low[i], lim_high[j],pre))
    #recall
    re = c(cM[1,1]/colSums(cM)[1], cM[2,2]/colSums(cM)[2], cM[3,3]/colSums(cM)[3])
    recall_all <- rbind(recall_all, c(lim_low[i], lim_high[j],lim_low[i], lim_high[j],re))
  }
}
colnames(recall_all) <- c("low_atac","high_atac", "low_wgs", "high_wgs", "recall_loss","recall_normal", "recall_gain")
colnames(precision_all) <- c("low_atac","high_atac", "low_wgs", "high_wgs", "precision_loss","precision_normal", "precision_gain")

#F1 scores per class:
F1_loss <- 2*precision_all[,5] * recall_all[,5] / (precision_all[,5] + recall_all[,5])
F1_normal <- 2*precision_all[,6] * recall_all[,6] / (precision_all[,6] + recall_all[,6])
F1_gain <- 2*precision_all[,7] * recall_all[,7] / (precision_all[,7] + recall_all[,7])
all<-cbind(precision_all,recall_all[,5:7],F1_loss,F1_normal,F1_gain)



#############################################################################################
#############################################################################################
#############################################################################################
#### plots:
#here come the limits used for the plot:
  #limits[1]=low threshold to call losses; 
  #limits[2]=high threshold to call gains
  #in between low and high it is the normal state
limits = c(1.8,2.1)
path1="path-to-plots/"

#PLOT:
#discretize the atac and wgs results into loss, normal, and gain given the low and high thresholds:
atac2 <- seg$scATAC
atac2[which(as.numeric(atac2) < as.numeric(limits[1]))] <- "loss"
atac2[which(as.numeric(atac2) > as.numeric(limits[2]))] <- "gain"
atac2[which(as.numeric(atac2) > as.numeric(limits[1]) & as.numeric(atac2) < as.numeric(limits[2]))] <- "normal"
wgs2 <- seg$scDNA
wgs2[which(as.numeric(wgs2) < as.numeric(limits[1]))] <- "loss"
wgs2[which(as.numeric(wgs2) > as.numeric(limits[2]))] <- "gain"
wgs2[which(as.numeric(wgs2) > as.numeric(limits[1]) & as.numeric(wgs2) < as.numeric(limits[2]))] <- "normal"
#attribute the loss, normal and gain to numerical values to plot them next to the signal:
wgs2[which(wgs2=="normal")] <- 2.3
wgs2[which(wgs2=="gain")] <- 3.3
wgs2[which(wgs2=="loss")] <- 1.3
atac2[which(atac2=="normal")] <- 2.0
atac2[which(atac2=="gain")] <- 3.0
atac2[which(atac2=="loss")] <- 1.0

#### plot all together
pdf(paste0(path1,"GenomePlots.pdf"), height=6, width=6)
layout(matrix(c(1,2,3),3,1,byrow=TRUE), heights=c(5,5,2.5))
par(mar = c(4, 4, 0.1, 0.1)) 

#plot the wgs profile genome-wide as a line
plot(seg$scDNA, type = 'l',col="cyan",ylim=c(1,3), ylab="wgs", xlab="")
abline(h=limits[1],col="lightblue")#low limit
abline(h=limits[2],col="lightblue")#high limit
#plot the atac profile genome-wide as a line
plot(seg$scATAC,type="l",col="orange", ylab="atac", xlab="")
abline(h=limits[1],col="pink")
abline(h=limits[2],col="pink")
#plot the segments of gain and loss for each modality:
plot(wgs2,col="darkcyan",ylim=c(0.8,3.7),yaxt="n", xlab="bin index", ylab="somy", cex=.4)
axis(2, at=c(1, 2, 3))
points(atac2,col="darkorange", cex=.4)
dev.off()



