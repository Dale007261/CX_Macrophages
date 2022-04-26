####Differential expression testing
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse, lib.loc = "D:/R-4.0.3/library")
library(ggrepel, lib.loc = "D:/R-4.0.3/library")
library("DESeq2")
library(EnhancedVolcano)
library(dplyr)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
#lpsvscon----
options(stringsAsFactors = FALSE) #禁止chr转成factor
Countlpscon<-read.csv("D:/CX_Macrophages/input/LPSvsCon.csv",header = T,sep = ",")
data2 <- Countlpscon[!duplicated(Countlpscon$gene_name),] 
rownames(data2)<-data2$gene_name
data2$gene_name<-NULL
group <- factor(c(2,2,2,1,1,1))
coldata <- data.frame(row.names = colnames(data2), group)
#PCA
data3<-t(data2)
pca.results <- prcomp(data3, center = TRUE, scale. = FALSE)
#colors
mycol <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
library(ggplot2)
library(plyr)
library(ggord)
library(yyplot)
ggord(pca.results, grp_in = coldata$group, repel=TRUE,
      ellipse = FALSE, 
      size = 2,
      alpha=0.5,
      cols = mycol[1:length(unique(coldata$group))],
      arrow = NULL,txt = NULL) + 
  theme(panel.grid =element_blank()) + 
  geom_ord_ellipse(ellipse_pro = .95,
                   size=1.5,
                   lty=1 )
dds <- DESeqDataSetFromMatrix(countData = data2,
                              colData = coldata,
                              design= ~ group )
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)
zzzzlps<-as.data.frame(res)
#IFNvsCon----
IFNvsCon<-read.csv("D:/CX_Macrophages/input/IFNvsCon.csv",header = T,sep = ",")
data2 <- IFNvsCon[!duplicated(IFNvsCon$gene_name),] 
rownames(data2)<-data2$gene_name
data2$gene_name<-NULL
group <- factor(c(2,2,2,1,1,1))
coldata <- data.frame(row.names = colnames(data2), group)

dds <- DESeqDataSetFromMatrix(countData = data2,
                              colData = coldata,
                              design= ~ group )
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)
zzzzifn<-as.data.frame(res)
#IFNCXvsIFN----
IFNCXvsIFN<-read.csv("D:/CX_Macrophages/input/IFNCXvsIFN.csv",header = T,sep = ",")
data2 <- IFNCXvsIFN[!duplicated(IFNCXvsIFN$gene_name),] 
rownames(data2)<-data2$gene_name
data2$gene_name<-NULL
group <- factor(c(1,1,1,2,2,2))
coldata <- data.frame(row.names = colnames(data2), group)

dds <- DESeqDataSetFromMatrix(countData = data2,
                              colData = coldata,
                              design= ~ group )
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)
zzzz<-as.data.frame(res)
#LPSCXvsLPS----
LPSCXvsLPS<-read.csv("D:/CX_Macrophages/input/LPSCXvsLPS.csv",header = T,sep = ",")
data2 <- LPSCXvsLPS[!duplicated(LPSCXvsLPS$gene_name),] 
rownames(data2)<-data2$gene_name
data2$gene_name<-NULL
group <- factor(c(2,2,2,1,1,1))
coldata <- data.frame(row.names = colnames(data2), group)

dds <- DESeqDataSetFromMatrix(countData = data2,
                              colData = coldata,
                              design= ~ group )
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)
zzzz<-as.data.frame(res)