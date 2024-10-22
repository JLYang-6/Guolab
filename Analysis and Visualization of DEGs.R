rt=read.table("Expression matrix file.txt",header = F,sep = "\t",check.names = F,quote = "")
rt=t(rt)
rt=rt[order(rt[,1]),]
rt=t(rt)
rt=rt[-1,-c(42:76)]
row.names(rt)=rt[,1]
colnames(rt)=rt[1,]
rt=rt[-1,-1]
rt1=matrix(as.numeric(rt),nrow=nrow(rt))
row.names(rt1)=row.names(rt)
colnames(rt1)=colnames(rt)
rt=rt1

library(limma)
rt=normalizeBetweenArrays(as.matrix(rt),method="scale")

class <- c(rep("CON",40),rep("PE",19))
design <- model.matrix(~factor(class))
colnames(design) <- c("CON","PE")
fit <- lmFit(rt,design)
fit2 <- eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=200000)

diffLab <- allDiff[with(allDiff, ((logFC> 0.5 | logFC< (-0.5)) & adj.P.Val < 0.05 )), ]
diffExpLevel <- rt[rownames(diffLab),]
qvalue=allDiff[rownames(diffLab),]$adj.P.Val
diffExpQvalue=cbind(qvalue,diffExpLevel)

library(pheatmap)
conNum=40                      
treatNum=19                      
geneFile="intersectGenes£¨top genes£©.txt"    
hmExp=diffExpLevel[as.vector(geneRT[,1]),]
hmExp=log2(hmExp+0.1)
Type=c(rep("Con",conNum),rep("PE",treatNum))
names(Type)=colnames(hmExp)
Type=as.data.frame(Type)
ann_colors <- list(Type = c(Con = '#1B19197F',PE = '#6318797F'))
pdf()
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue2", "white", "red2"))(50),
         annotation_colors=ann_colors,
         border_color = "black",
         cluster_cols =F,
         cluster_rows = F,
         show_colnames = F,
         show_rownames = T,
         scale="row",
         fontsize = 12,
         fontsize_row=8,
         fontsize_col=10)
dev.off()

library(ggplot2) 
library(ggrepel)
library(Seurat) 
logFCfilter=0.5
fdrFilter=0.01
Significant=ifelse((allDiff$adj.P.Val<fdrFilter & abs(allDiff$logFC)>logFCfilter), ifelse(allDiff$logFC>logFCfilter,"Up in preeclampsia","Down in preeclampsia"), "Insignificant")
p = ggplot(allDiff, aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(shape=21,aes(col=Significant ))+
  scale_color_manual("type",values=c("Down in preeclampsia"="#3333CC","Insignificant"="#CCCCCC","Up in preeclampsia"="#CC0000"))+
  labs(title = "")+
  scale_x_continuous("log2 (FoldChange)",limits=c(-2.5,2.5))+
  scale_y_continuous("-log10 (Adjusted p value)",limits = c(0,7.5))+
  geom_vline(xintercept=c(-0.5,0.5),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.01),lty=2,col="black",lwd=0.5)+
  theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))
p+theme_bw()
plot1 <- LabelPoints(plot = p, points = "MXRA5", label= "MXRA5",repel = FALSE,xnudge = -0.2,ynudge = 0.2)
dev.off()