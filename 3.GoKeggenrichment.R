####Go_Kegg_enrichment
#LPSvsCon
LPSvsConDEG<-read.csv("D:/CX_Macrophages/output/LPSvsConDEG.csv",header = T,sep = ",")
gene<-LPSvsConDEG[,1]
ego <- enrichGO(gene         = object,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1)
write.csv(ego@result,"D:/CX_Macrophages/output/result/GO_LPSvsCon.csv")
gene <- as.character(gene)
object1 = bitr(gene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
kk <- enrichKEGG(gene         = object1$ENTREZID,
                 organism     = 'mmu',
                 pvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
kk<-kk@result
write.csv(kk,"D:/CX_Macrophages/output/result/kegg_LPSvsCon.csv")
#IFNvsCon
IFNvsConDEG<-read.csv("D:/CX_Macrophages/output/IFNvsConDEG.csv",header = T,sep = ",")
gene<-LPSvsConDEG[,1]
ego <- enrichGO(gene         = object,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1)
write.csv(ego@result,"D:/CX_Macrophages/output/result/GO_IFNvsCon.csv")
gene <- as.character(gene)
object1 = bitr(gene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
kk <- enrichKEGG(gene         = object1$ENTREZID,
                 organism     = 'mmu',
                 pvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
kk<-kk@result
write.csv(kk,"D:/CX_Macrophages/output/result/kegg_IFNvsCon.csv")
#IFNvsIFNCX
IFNvsIFNCXDEG<-read.csv("D:/CX_Macrophages/output/IFNvsIFNCXDEG.csv",header = T,sep = ",")
gene<-IFNvsIFNCXDEG[,1]
ego <- enrichGO(gene         = object,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1)
write.csv(ego@result,"D:/CX_Macrophages/output/result/GO_IFNvsIFNCX.csv")
object <- as.character(object)
object1 = bitr(object, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
kk <- enrichKEGG(gene         = object1$ENTREZID,
                 organism     = 'mmu',
                 pvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
kk<-kk@result
write.csv(kk,"D:/CX_Macrophages/output/result/kegg_IFNvsIFNCX.csv")
#lpsvslpscx
LPSCXvsLPSDEG<-read.csv("D:/CX_Macrophages/output/LPSCXvsLPSDEG.csv",header = T,sep = ",")
gene<-LPSCXvsLPSDEG[,1]
ego <- enrichGO(gene         = object,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1)
write.csv(ego@result,"D:/CX_Macrophages/output/result/GO_LPSvsLPSCX.csv")
gene <- as.character(gene)
object1 = bitr(gene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
kk <- enrichKEGG(gene         = object1$ENTREZID,
                 organism     = 'mmu',
                 pvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
kk<-kk@result
write.csv(kk,"D:/CX_Macrophages/output/result/kegg_LPSvsLPSCX.csv")