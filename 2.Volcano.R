####Volcano
library(EnhancedVolcano)
x<-zzzzlps
gene<-c("Il1a","Il1b","Il6","Nos2","Ccl2","Ccl3","Cxcl2","Ccl12")
p1<-EnhancedVolcano(x,
                    lab = rownames(x),
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    selectLab = c("Il1a","Il1b","Il6","Nos2","Ccl2","Ccl3","Cxcl2","Ccl12"),
                    xlab = bquote(~Log[2]~ 'fold change'),
                    pCutoff = 0.05,
                    FCcutoff = 1.0,
                    pointSize = 2,
                    labSize = 4.0,
                    labCol = 'black',
                    labFace = 'bold',
                    boxedLabels = F,
                    colAlpha = 4/5,
                    legendPosition = 'right',
                    legendLabSize = 14,
                    legendIconSize = 4.0,
                    drawConnectors = T,
                    widthConnectors = 0.05,
                    colConnectors = 'black')
p1
x<-zzzzifn
p2<-EnhancedVolcano(x,
                    lab = rownames(x),
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    selectLab = c("Ccl4","Ccl6","Slamf1"),
                    xlab = bquote(~Log[2]~ 'fold change'),
                    pCutoff = 0.05,
                    FCcutoff = 1.0,
                    pointSize = 2.0,
                    labSize = 4.0,
                    labCol = 'black',
                    labFace = 'bold',
                    boxedLabels = F,
                    colAlpha = 4/5,
                    legendPosition = 'right',
                    legendLabSize = 14,
                    legendIconSize = 4.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.05,
                    colConnectors = 'black')
p2
