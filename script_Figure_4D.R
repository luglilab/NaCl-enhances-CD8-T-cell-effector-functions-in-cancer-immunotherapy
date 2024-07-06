###################
###################
# Code Figure: 4D
###################
###################

## FIGURE 4D
### Heatmap bulk rna glutamine ###


module load conda/anaconda3
source activate bulkRNAseq

R
set.seed(123)

library(edgeR)

# load raw matrix
counts <- read.table(file="Dataset_3_raw_count_matrix.txt", header= TRUE, stringsAsFactors=FALSE, sep="\t", row.names = 1)

# load gene annotation
anno <- read.table(file= "gene_annotation.txt", header= TRUE, stringsAsFactors=FALSE, sep="\t")

# order samples
counts <- counts[c("Lib_CSC-37", "Lib_CSC-43", "Lib_CSC-49", "Lib_CSC-55", "Lib_CSC-39", "Lib_CSC-44", "Lib_CSC-51", "Lib_CSC-57", "Lib_CSC-41", "Lib_CSC-47", "Lib_CSC-53", "Lib_CSC-59", "Lib_CSC-42", "Lib_CSC-48", "Lib_CSC-54", "Lib_CSC-60")]

Treat <- factor(c(rep("CT",4), rep("CT_no_glutamine",4), rep("NaCl",4), rep("NaCl_no_glutamine",4)))

Sample <- factor(rep(c("A", "B", "C", "D"),4))

dge_data <- DGEList(counts, group=Treat)

keep <- filterByExpr(dge_data, group=Treat)

dge_data <- dge_data[keep,,keep.lib.sizes=FALSE]

dge_data <- calcNormFactors(dge_data)
dge_data$samples

design <- model.matrix(~Sample+Treat)
dge_data <- estimateDisp(dge_data,design,robust=TRUE)

matrice_cpm <- cpm(dge_data)
matrice_log2cpm <- log2(matrice_cpm+1)

matrice_log2cpm_anno <- merge(anno,matrice_log2cpm,by.x="gene_ID",by.y="row.names",sort=FALSE)

# remove extra anno
matrice_log2cpm_anno_v2 <- matrice_log2cpm_anno[c("gene_Symbol","Lib_CSC-37", "Lib_CSC-43", "Lib_CSC-49", "Lib_CSC-55", "Lib_CSC-39", "Lib_CSC-44", "Lib_CSC-51", "Lib_CSC-57", "Lib_CSC-41", "Lib_CSC-47", "Lib_CSC-53", "Lib_CSC-59", "Lib_CSC-42", "Lib_CSC-48", "Lib_CSC-54", "Lib_CSC-60")]

matrice_log2cpm_anno_v2 <- matrice_log2cpm_anno_v2 %>% distinct(gene_Symbol, .keep_all = TRUE)
rownames(matrice_log2cpm_anno_v2) <- matrice_log2cpm_anno_v2 $gene_Symbol
matrice_log2cpm_anno_v2 <- matrice_log2cpm_anno_v2[,-1]

experiment <- as.data.frame(cbind(c("Lib_CSC-37", "Lib_CSC-43", "Lib_CSC-49", "Lib_CSC-55", "Lib_CSC-39", "Lib_CSC-44", "Lib_CSC-51", "Lib_CSC-57", "Lib_CSC-41", "Lib_CSC-47", "Lib_CSC-53", "Lib_CSC-59", "Lib_CSC-42", "Lib_CSC-48", "Lib_CSC-54", "Lib_CSC-60"),c("Lib_CSC-37", "Lib_CSC-43", "Lib_CSC-49", "Lib_CSC-55", "Lib_CSC-39", "Lib_CSC-44", "Lib_CSC-51", "Lib_CSC-57", "Lib_CSC-41", "Lib_CSC-47", "Lib_CSC-53", "Lib_CSC-59", "Lib_CSC-42", "Lib_CSC-48", "Lib_CSC-54", "Lib_CSC-60"),c(rep("CT",4), rep("CT_no_glutamine",4), rep("NaCl",4), rep("NaCl_no_glutamine",4))))
colnames(experiment) <- c("Samples", "Sample", "Condition")

class <- as.factor(experiment$Condition) 
identical(colnames(matrice_log2cpm_anno_v2), as.character(experiment$Sample))

matrx <- NULL                                                        
for (i in 1:length(rownames(matrice_log2cpm_anno_v2)))                   
{vector<-anova(lm(as.numeric(matrice_log2cpm_anno_v2[i,]) ~ class))      
 vett <- vector$"Pr(>F)"                                             
 matrx <- c(matrx,vett[1])}                                          
 matrx2 <- as.matrix(matrx)                                          
 namesMatrice <- rownames(matrice_log2cpm_anno_v2)                          
 rownames(matrx2) <- namesMatrice       

DEGs <- as.matrix(matrx2[matrx2<0.05,])

require(multtest)
procs <- c("BH","Bonferroni")
res <- mt.rawp2adjp(matrx2, procs)
adjp <- res$adjp[order(res$index), ]


rownames(adjp) <- rownames(matrx2)
write.table(adjp, "anova_results.txt", sep="\t")

anova_test <- read.table(file="anova_results.txt", header= TRUE, stringsAsFactors=FALSE, sep="\t", row.names = 1)

# filter Bonferroni < 0.05
gene_significant <- anova_test[anova_test$Bonferroni < 0.05,]
gene <- rownames(gene_significant)

matrice_log2cpm_anno_v3 <- matrice_log2cpm_anno_v2[gene,]
write.table(matrice_log2cpm_anno_v3, "Raw_data_FIG4_D.txt", sep="\t", col.name = TRUE, quote = FALSE)

library(gplots)
matrice_heat <- matrice_log2cpm_anno_v3

heat.colors <- colorRampPalette(c("violet", "black", "gold"))(256)
scale.data <- as.matrix((matrice_heat-apply(matrice_heat,1,mean))/apply(matrice_heat,1,sd))

# FIGURE 4D

pdf("FIG4_D.pdf", height=40, width=22)
heatmap.2(scale.data,Colv = T,
            dendrogram="both",scale="none",colsep=c(0),
            sepcolor="white",sepwidth=c(0.01),margins = c(13,13),
            col=heat.colors,trace="none",labRow = rownames(matrice_heat),labCol = colnames(matrice_heat),
            key=T,keysize=0.8,density.info="none",symkey=TRUE, key.xlab=c("log2ratio"),cexRow=0.45,ColSideColors=c(rep("blue",4),rep("red",4),rep("lightpink1",4),rep("lightblue1",4)))
            legend("topright",legend=c("CT","CT_no_glut","NaCl","NaCl_no_glut"),col=c("blue", "red", "lightpink1", "lightblue1"),pch=c(15,15,15,15),cex=4,pt.cex=5,bty = "n")
par(cex.main=0.5)
dev.off()