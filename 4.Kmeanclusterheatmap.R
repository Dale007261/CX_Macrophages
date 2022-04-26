library(ggplot2) 		#heatmap plotting
library(tidyverse)  		# data manipulation
library(cluster)    		# clustering algorithms
setwd("D:/CX_Macrophages/input")
data<-read.csv("LPSIFNFPKM.csv",header = T)
data <- data[!duplicated(data$gene_name),] 
rownames(data)<-data$gene_name
data$gene_name<-NULL
y <- as.matrix(data)
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- apply(y,2,fpkmToTpm)
tpms[1:3,]
colSums(tpms)

boxplot(y,las=2)
exprSet <- tpms
group_list=c(rep('control',3),rep('IFN',3),rep('LPS',3))
group_list <- factor(group_list,levels = c("control","IFN","LPS"),ordered = F)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
y<-exprSet
IFNDEG<-read.csv("D:/CX_Macrophages/output/IFNvsConDEG.csv",header=T,row.names = 1)
IFNDEG<-rownames(IFNDEG)
LPSvsConDEG<-read.csv("D:/CX_Macrophages/output/LPSvsConDEG.csv",header=T,row.names = 1)
LPSDEG<-rownames(LPSvsConDEG)
Alldeg<-c(IFNDEG,LPSDEG)
Alldeg<-unique(Alldeg)
y<-as.data.frame(y)
y<-filter(y, row.names(y) %in% Alldeg)
km <- kmeans(t(scale(t(y))), 2)
kmcluster<-km[["cluster"]]
hr <- hclust(as.dist(1 - cor(t(y), method = "pearson")), method = "complete")
hc <- hclust(as.dist(1 - cor(y, method = "spearman")), method = "complete")
mycl <-  cutree(hr, h = max(hr$height) / 1.01)
mycolhc <-   rainbow(length(unique(mycl)), start = 0.1, end = 0.9)
mycolhc <- mycolhc[as.vector(mycl)]
output<-cbind(y,kmcluster)
mycol <-  colorpanel(75, "blue", "white", "yellow") # or try redgreen(75)
write.csv(output,"Kmeancluster.csv")
y<-as.matrix(y)
library(gplots)
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol,
          scale="row", density.info="none", trace="none", 
          RowSideColors=as.character(kmcluster))
#IFNCXLPSCX
library(ggplot2) 		#heatmap plotting
library(tidyverse)  		# data manipulation
library(cluster)    		# clustering algorithms
setwd("D:/CX_Macrophages/input")
data<-read.csv("LPSIFNCXFPKM.csv",header = T)
data <- data[!duplicated(data$gene_name),] 
rownames(data)<-data$gene_name
data$gene_name<-NULL
y <- as.matrix(data)
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- apply(y,2,fpkmToTpm)
tpms[1:3,]
colSums(tpms)

boxplot(y,las=2)
exprSet <- tpms
group_list=c(rep('IFNcx',3),rep('LPScx',3),rep('IFN',3),rep('LPS',3),)
group_list <- factor(group_list,levels = c("IFNcx","LPScx","IFN","LPS"),ordered = F)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
LPSCXvsLPSDEG<-read.csv("D:/CX_Macrophages/output/LPScx/LPSCXvsLPSDEG.csv",header=T,row.names = 1)
LPSCXvsLPSDEG<-rownames(LPSCXvsLPSDEG)
IFNCXvsIFNDEG<-read.csv("D:/CX_Macrophages/output/IFNcx/IFNCXvsIFNDEG.csv",header=T,row.names = 1)
IFNCXvsIFNDEG<-rownames(IFNCXvsIFNDEG)
Alldeg<-c(LPSCXvsLPSDEG,IFNCXvsIFNDEG)
Alldeg<-unique(Alldeg)
y<-as.data.frame(y)
y<-filter(y, row.names(y) %in% Alldeg)
km <- kmeans(t(scale(t(y))), 2)
kmcluster<-km[["cluster"]]
hr <- hclust(as.dist(1 - cor(t(y), method = "pearson")), method = "complete")
hc <- hclust(as.dist(1 - cor(y, method = "spearman")), method = "complete")
mycl <-  cutree(hr, h = max(hr$height) / 1.01)
mycolhc <-   rainbow(length(unique(mycl)), start = 0.1, end = 0.9)
mycolhc <- mycolhc[as.vector(mycl)]
output<-cbind(y,kmcluster)
mycol <-  colorpanel(75, "blue", "white", "yellow") # or try redgreen(75)
write.csv(output,"Kmeancluster.csv")
y<-as.matrix(y)
library(gplots)
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol,
          scale="row", density.info="none", trace="none", 
          RowSideColors=as.character(kmcluster))
