###################
###################
# Code Figure: 3A, 3D, 3F
###################
###################

## FIGURE 3A
### Luoma Tumor ###

module load conda/anaconda3
source activate scRNAseq

R
set.seed(123)

# Signature from Tumor Luoma
gene_signature <- read.table("Luoma_Tumor_CD8_signature.txt", header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep='\t')
dim(gene_signature)
# 677   7

# Signature CD8 HSD/NSD
final_summary <- read.table("SC_Mouse_cluster4_HSDvsNSD_convertedHuman.txt", header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep='\t')
dim(final_summary)
# 243   1

###### pvalue ######
total_enrichment <- NULL
total_label <- NULL

 for (i in 1:ncol(final_summary)) {
    print(paste("enrichment for",i))
    anti_probesets <- final_summary[,i]
    anti_symbols <- anti_probesets[!is.na(anti_probesets)]
    anti_symbols<-anti_symbols[anti_symbols!=""]

    anti_symbols<-unique(anti_symbols)

    ### enrichment test
    pvalues_symbols=NULL

    for (j in 1:dim(gene_signature)[2]) {

      sig_symbols<-gene_signature[,j]
      sig_symbols<-sig_symbols[sig_symbols!=""]
      sig_symbols<-unique(sig_symbols)

      ### reduce to background
      sig_symbols<-sig_symbols


      ### ipergeometrico sui symbols
      common_symbols<- sig_symbols[sig_symbols %in% anti_symbols]
      x=length(common_symbols)
      if (x>0) {
      m=length(sig_symbols)               #white balls in the urn.
      n=16566 - m      #black balls in the urn
      k=length(anti_symbols)                    #balls drawn from the urn
      pvalue<-phyper( x-1,m,n,k,lower.tail=F)
      }  else {pvalue=NA}
      pvalues_symbols=c(pvalues_symbols,pvalue)
    } # for signatures

    total_enrichment<-cbind(total_enrichment, pvalues_symbols)
    total_label<-c(total_label,colnames(final_summary)[i])

  }

  colnames(total_enrichment)<-total_label
rownames(total_enrichment)<-colnames(gene_signature)

pval_table <- total_enrichment

###### core_gene ######
total_enrichment <- NULL
total_label <- NULL

 for (i in 1:ncol(final_summary)) {
    print(paste("enrichment for",i))
    anti_probesets <- final_summary[,i]
    anti_symbols <- anti_probesets[!is.na(anti_probesets)]
    anti_symbols<-anti_symbols[anti_symbols!=""]

    anti_symbols<-unique(anti_symbols)

    ### enrichment test

    for (j in 1:dim(gene_signature)[2]) {

      sig_symbols<-gene_signature[,j]
      sig_symbols<-sig_symbols[sig_symbols!=""]
      sig_symbols<-unique(sig_symbols)

      ### reduce to background
      sig_symbols<-sig_symbols


      ### ipergeometrico sui symbols
      common_symbols<- sig_symbols[sig_symbols %in% anti_symbols]
      core_gene=length(common_symbols)
          total_enrichment<-cbind(total_enrichment, core_gene)
    } # for signatures

    total_label<-c(total_label,colnames(final_summary)[i])

  }


core_gene_table <- total_enrichment

pval_table <- as.data.frame(pval_table)
pval_shape <- data.frame(rows = row.names(pval_table), stack(pval_table))
core_gene <- as.data.frame(t(core_gene_table))

final_baloon_table <- cbind(pval_shape$rows,core_gene$V1,pval_shape$values,as.character(pval_shape$ind))
colnames(final_baloon_table) <- c("signature","core_gene","pval","group")
write.table(final_baloon_table, "raw_data_FIG3_A.txt",sep="\t", quote = FALSE, row.names = FALSE)

# reload table to plot correctly
library(ggplot2)
library(ggpubr)
data <- read.table("raw_data_FIG3_A.txt", sep="\t", header=T)
data <- data[8:14,1:4]
data <- as.data.frame(data)
data$pval <- -log10(data$pval)

data$pval[data$pval < 1.30103 ] <- NA

col_scale <- colorRampPalette(c("blue","pink", "red", "darkred"))(10)


myplot <- ggballoonplot(data,x = "group", y = "signature", size = "core_gene", fill = "pval") + scale_fill_gradientn(colors = col_scale, na.value="white")  + labs(size="Core gene",fill="-log10(pval)") + scale_size_continuous(limits=c(0,max(data$core_gene)))

# FIGURE 3A
pdf("FIG3_A.pdf")
myplot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


####################


## FIGURE 3D
### Luoma Blood ###

module load conda/anaconda3
source activate scRNAseq

R
set.seed(123)

# Signature from Blood Luoma
gene_signature <- read.table("Luoma_Blood_CD8_signature.txt", header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep='\t')
dim(gene_signature)
# 535   9

# Signature CD8 HSD/NSD
final_summary <- read.table("SC_Mouse_cluster4_HSDvsNSD_convertedHuman.txt", header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep='\t')
dim(final_summary)
# 243   1

###### pvalue ######
total_enrichment <- NULL
total_label <- NULL

 for (i in 1:ncol(final_summary)) {
    print(paste("enrichment for",i))
    anti_probesets <- final_summary[,i]
    anti_symbols <- anti_probesets[!is.na(anti_probesets)]
    anti_symbols<-anti_symbols[anti_symbols!=""]

    anti_symbols<-unique(anti_symbols)

    ### enrichment test
    pvalues_symbols=NULL

    for (j in 1:dim(gene_signature)[2]) {

      sig_symbols<-gene_signature[,j]
      sig_symbols<-sig_symbols[sig_symbols!=""]
      sig_symbols<-unique(sig_symbols)

      ### reduce to background
      sig_symbols<-sig_symbols


      ### ipergeometrico sui symbols
      common_symbols<- sig_symbols[sig_symbols %in% anti_symbols]
      x=length(common_symbols)
      if (x>0) {
      m=length(sig_symbols)               #white balls in the urn.
      n=16566 - m      #black balls in the urn
      k=length(anti_symbols)                    #balls drawn from the urn
      pvalue<-phyper( x-1,m,n,k,lower.tail=F)
      }  else {pvalue=NA}
      pvalues_symbols=c(pvalues_symbols,pvalue)
    } # for signatures

    total_enrichment<-cbind(total_enrichment, pvalues_symbols)
    total_label<-c(total_label,colnames(final_summary)[i])

  }

  colnames(total_enrichment)<-total_label
rownames(total_enrichment)<-colnames(gene_signature)

pval_table <- total_enrichment

###### core_gene ######
total_enrichment <- NULL
total_label <- NULL

 for (i in 1:ncol(final_summary)) {
    print(paste("enrichment for",i))
    anti_probesets <- final_summary[,i]
    anti_symbols <- anti_probesets[!is.na(anti_probesets)]
    anti_symbols<-anti_symbols[anti_symbols!=""]

    anti_symbols<-unique(anti_symbols)

    ### enrichment test

    for (j in 1:dim(gene_signature)[2]) {

      sig_symbols<-gene_signature[,j]
      sig_symbols<-sig_symbols[sig_symbols!=""]
      sig_symbols<-unique(sig_symbols)

      ### reduce to background
      sig_symbols<-sig_symbols


      ### ipergeometrico sui symbols
      common_symbols<- sig_symbols[sig_symbols %in% anti_symbols]
      core_gene=length(common_symbols)
          total_enrichment<-cbind(total_enrichment, core_gene)
    } # for signatures

    total_label<-c(total_label,colnames(final_summary)[i])

  }


core_gene_table <- total_enrichment

pval_table <- as.data.frame(pval_table)
pval_shape <- data.frame(rows = row.names(pval_table), stack(pval_table))
core_gene <- as.data.frame(t(core_gene_table))

final_baloon_table <- cbind(pval_shape$rows,core_gene$V1,pval_shape$values,as.character(pval_shape$ind))
colnames(final_baloon_table) <- c("signature","core_gene","pval","group")
write.table(final_baloon_table, "raw_data_FIG3_D.txt",sep="\t", quote = FALSE, row.names = FALSE)

# reload table to plot correctly
library(ggplot2)
library(ggpubr)
data <- read.table("raw_data_FIG3_D.txt", sep="\t", header=T)
data <- data[,1:4]
data <- as.data.frame(data)
data$pval <- -log10(data$pval)

data$pval[data$pval < 1.30103 ] <- NA

col_scale <- colorRampPalette(c("blue","pink", "red", "darkred"))(10)


myplot <- ggballoonplot(data,x = "group", y = "signature", size = "core_gene", fill = "pval") + scale_fill_gradientn(colors = col_scale, na.value="white")  + labs(size="Core gene",fill="-log10(pval)") + scale_size_continuous(limits=c(0,max(data$core_gene)))

# FIGURE 3D
pdf("FIG3_D.pdf")
myplot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


####################


## FIGURE 3F TOP
### Fairfax Response ###

module load conda/anaconda3
source activate scRNAseq

R
set.seed(123)

# Signature from Fairfax
gene_signature <- read.table("Fairfax_signature.txt", header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep='\t')
dim(gene_signature)
# 2384   2

# Signature CD8 HSD/NSD
final_summary <- read.table("SC_Mouse_cluster4_HSDvsNSD_convertedHuman.txt", header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep='\t')
dim(final_summary)
# 243   1

###### pvalue ######
total_enrichment <- NULL
total_label <- NULL

 for (i in 1:ncol(final_summary)) {
    print(paste("enrichment for",i))
    anti_probesets <- final_summary[,i]
    anti_symbols <- anti_probesets[!is.na(anti_probesets)]
    anti_symbols<-anti_symbols[anti_symbols!=""]

    anti_symbols<-unique(anti_symbols)

    ### enrichment test
    pvalues_symbols=NULL

    for (j in 1:dim(gene_signature)[2]) {

      sig_symbols<-gene_signature[,j]
      sig_symbols<-sig_symbols[sig_symbols!=""]
      sig_symbols<-unique(sig_symbols)

      ### reduce to background
      sig_symbols<-sig_symbols


      ### ipergeometrico sui symbols
      common_symbols<- sig_symbols[sig_symbols %in% anti_symbols]
      x=length(common_symbols)
      if (x>0) {
      m=length(sig_symbols)               #white balls in the urn.
      n=16566 - m      #black balls in the urn
      k=length(anti_symbols)                    #balls drawn from the urn
      pvalue<-phyper( x-1,m,n,k,lower.tail=F)
      }  else {pvalue=NA}
      pvalues_symbols=c(pvalues_symbols,pvalue)
    } # for signatures

    total_enrichment<-cbind(total_enrichment, pvalues_symbols)
    total_label<-c(total_label,colnames(final_summary)[i])

  }

  colnames(total_enrichment)<-total_label
rownames(total_enrichment)<-colnames(gene_signature)

pval_table <- total_enrichment

###### core_gene ######
total_enrichment <- NULL
total_label <- NULL

 for (i in 1:ncol(final_summary)) {
    print(paste("enrichment for",i))
    anti_probesets <- final_summary[,i]
    anti_symbols <- anti_probesets[!is.na(anti_probesets)]
    anti_symbols<-anti_symbols[anti_symbols!=""]

    anti_symbols<-unique(anti_symbols)

    ### enrichment test

    for (j in 1:dim(gene_signature)[2]) {

      sig_symbols<-gene_signature[,j]
      sig_symbols<-sig_symbols[sig_symbols!=""]
      sig_symbols<-unique(sig_symbols)

      ### reduce to background
      sig_symbols<-sig_symbols


      ### ipergeometrico sui symbols
      common_symbols<- sig_symbols[sig_symbols %in% anti_symbols]
      core_gene=length(common_symbols)
          total_enrichment<-cbind(total_enrichment, core_gene)
    } # for signatures

    total_label<-c(total_label,colnames(final_summary)[i])

  }


core_gene_table <- total_enrichment

pval_table <- as.data.frame(pval_table)
pval_shape <- data.frame(rows = row.names(pval_table), stack(pval_table))
core_gene <- as.data.frame(t(core_gene_table))

final_baloon_table <- cbind(pval_shape$rows,core_gene$V1,pval_shape$values,as.character(pval_shape$ind))
colnames(final_baloon_table) <- c("signature","core_gene","pval","group")
write.table(final_baloon_table, "raw_data_FIG3_F_TOP.txt",sep="\t", quote = FALSE, row.names = FALSE)
Fairfax_result <- final_baloon_table

## FIGURE 3F BOTTOM
### Caushi Response ###

module load conda/anaconda3
source activate scRNAseq

R
set.seed(123)

# Signature from Caushi
gene_signature <- read.table("Caushi_CD8_tumor_MPR_Vs_nonMPR_signature.txt", header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep='\t')
dim(gene_signature)
# 85   2

# Signature CD8 HSD/NSD
final_summary <- read.table("SC_Mouse_cluster4_HSDvsNSD_convertedHuman.txt", header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep='\t')
dim(final_summary)
# 243   1

###### pvalue ######
total_enrichment <- NULL
total_label <- NULL

 for (i in 1:ncol(final_summary)) {
    print(paste("enrichment for",i))
    anti_probesets <- final_summary[,i]
    anti_symbols <- anti_probesets[!is.na(anti_probesets)]
    anti_symbols<-anti_symbols[anti_symbols!=""]

    anti_symbols<-unique(anti_symbols)

    ### enrichment test
    pvalues_symbols=NULL

    for (j in 1:dim(gene_signature)[2]) {

      sig_symbols<-gene_signature[,j]
      sig_symbols<-sig_symbols[sig_symbols!=""]
      sig_symbols<-unique(sig_symbols)

      ### reduce to background
      sig_symbols<-sig_symbols


      ### ipergeometrico sui symbols
      common_symbols<- sig_symbols[sig_symbols %in% anti_symbols]
      x=length(common_symbols)
      if (x>0) {
      m=length(sig_symbols)               #white balls in the urn.
      n=16566 - m      #black balls in the urn
      k=length(anti_symbols)                    #balls drawn from the urn
      pvalue<-phyper( x-1,m,n,k,lower.tail=F)
      }  else {pvalue=NA}
      pvalues_symbols=c(pvalues_symbols,pvalue)
    } # for signatures

    total_enrichment<-cbind(total_enrichment, pvalues_symbols)
    total_label<-c(total_label,colnames(final_summary)[i])

  }

  colnames(total_enrichment)<-total_label
rownames(total_enrichment)<-colnames(gene_signature)

pval_table <- total_enrichment

###### core_gene ######
total_enrichment <- NULL
total_label <- NULL

 for (i in 1:ncol(final_summary)) {
    print(paste("enrichment for",i))
    anti_probesets <- final_summary[,i]
    anti_symbols <- anti_probesets[!is.na(anti_probesets)]
    anti_symbols<-anti_symbols[anti_symbols!=""]

    anti_symbols<-unique(anti_symbols)

    ### enrichment test

    for (j in 1:dim(gene_signature)[2]) {

      sig_symbols<-gene_signature[,j]
      sig_symbols<-sig_symbols[sig_symbols!=""]
      sig_symbols<-unique(sig_symbols)

      ### reduce to background
      sig_symbols<-sig_symbols


      ### ipergeometrico sui symbols
      common_symbols<- sig_symbols[sig_symbols %in% anti_symbols]
      core_gene=length(common_symbols)
          total_enrichment<-cbind(total_enrichment, core_gene)
    } # for signatures

    total_label<-c(total_label,colnames(final_summary)[i])

  }


core_gene_table <- total_enrichment

pval_table <- as.data.frame(pval_table)
pval_shape <- data.frame(rows = row.names(pval_table), stack(pval_table))
core_gene <- as.data.frame(t(core_gene_table))

final_baloon_table <- cbind(pval_shape$rows,core_gene$V1,pval_shape$values,as.character(pval_shape$ind))
colnames(final_baloon_table) <- c("signature","core_gene","pval","group")
write.table(final_baloon_table, "raw_data_FIG3_F_BOTTOM.txt",sep="\t", quote = FALSE, row.names = FALSE)
Caushi_result <- final_baloon_table


Figure3_F <- rbind(Fairfax_result,Caushi_result)
write.table(Figure3_F, "FIG3_F.txt",sep="\t", quote = FALSE, row.names = FALSE)

