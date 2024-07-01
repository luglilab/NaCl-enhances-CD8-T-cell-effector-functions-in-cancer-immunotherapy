###################
###################
# Mouse SC-Rna Cd8 recluster
###################
###################

module load conda/anaconda3
source activate scRNAseq

R
set.seed(123)

library(Seurat)
library(dplyr)
library(qpcR)
library(ggplot2)

# start from first analysis data

matrix_raw <- read.table("raw_counts_marco.csv", sep="\t", header = TRUE, row.names = 1)
matrix_raw <- t(matrix_raw)
dim(matrix_raw)
# 16566 37525

metadata <- read.table("metadata_marco.csv", sep="\t", header=T, row.names=1)
dim(metadata)
# 37525    28

# add res 0.5 annotation

#0	Mono/Macro
#1	Macrophages
#2	Neutrophils
#3	NK cells
#4	CD8+ T cells
#5	Macrophages
#6	NK cells
#7	CD4+ T cells
#8	Neutrophils Arg1,2+
#9	Neutrophils
#10	Macrophages
#11	DCs
#12	Macrophages
#13	DCs

#    0    1    2    3    4    5    6    7    8    9   10   11   12   13 
# 7579 5017 5474 2909 2582 3162 1206 1120 1004 2200 1440  601 2908  323 

anno_desc_0.5 <- metadata$desc_0.5
metadata <- cbind(metadata,anno_desc_0.5)

metadata$anno_desc_0.5 <- gsub(x = metadata$anno_desc_0.5, pattern = "11", replacement = "DCs", fixed = TRUE)  
metadata$anno_desc_0.5 <- gsub(x = metadata$anno_desc_0.5, pattern = "10", replacement = "Macrophages", fixed = TRUE)  
metadata$anno_desc_0.5 <- gsub(x = metadata$anno_desc_0.5, pattern = "12", replacement = "Macrophages", fixed = TRUE)  
metadata$anno_desc_0.5 <- gsub(x = metadata$anno_desc_0.5, pattern = "13", replacement = "DCs", fixed = TRUE)  
metadata$anno_desc_0.5 <- gsub(x = metadata$anno_desc_0.5, pattern = "1", replacement = "Macrophages", fixed = TRUE)  
metadata$anno_desc_0.5 <- gsub(x = metadata$anno_desc_0.5, pattern = "0", replacement = "Mono_Macro", fixed = TRUE)  
metadata$anno_desc_0.5 <- gsub(x = metadata$anno_desc_0.5, pattern = "2", replacement = "Neutrophils", fixed = TRUE)  
metadata$anno_desc_0.5 <- gsub(x = metadata$anno_desc_0.5, pattern = "5", replacement = "Macrophages", fixed = TRUE)  
metadata$anno_desc_0.5 <- gsub(x = metadata$anno_desc_0.5, pattern = "9", replacement = "Neutrophils", fixed = TRUE)  
metadata$anno_desc_0.5 <- gsub(x = metadata$anno_desc_0.5, pattern = "3", replacement = "NK_cells", fixed = TRUE)  
metadata$anno_desc_0.5 <- gsub(x = metadata$anno_desc_0.5, pattern = "6", replacement = "NK_cells", fixed = TRUE)  
metadata$anno_desc_0.5 <- gsub(x = metadata$anno_desc_0.5, pattern = "8", replacement = "Neutrophils_arg", fixed = TRUE)  
metadata$anno_desc_0.5 <- gsub(x = metadata$anno_desc_0.5, pattern = "4", replacement = "CD8_Tcells", fixed = TRUE)  
metadata$anno_desc_0.5 <- gsub(x = metadata$anno_desc_0.5, pattern = "7", replacement = "CD4_Tcells", fixed = TRUE)  

#     CD4_Tcells      CD8_Tcells             DCs     Macrophages      Mono_Macro 
#           1120            2582             924           12527            7579 
#    Neutrophils Neutrophils_arg        NK_cells 
#           7674            1004            4115 

# Select non doublets cells
cellule_NOdoublets <- metadata[metadata$predicted_doublets == "False",]
dim(cellule_NOdoublets)
# 35932    28

matrix_raw_NOdoublets <- matrix_raw[,rownames(cellule_NOdoublets)]
dim(matrix_raw_NOdoublets)
# 16566 35932

data_cells <- CreateSeuratObject(counts = matrix_raw_NOdoublets, min.cells = 3, min.features = 200)
dim(data_cells)
# 16566 35932

data_cells <- AddMetaData(data_cells, metadata = cellule_NOdoublets)
head(data_cells@meta.data, 5)

saveRDS(data_cells, file = "SC_topo_NaCl_anno.rds")


##########

module load conda/anaconda3
source activate scRNAseq

R
set.seed(123)

library(Seurat)
library(dplyr)
library(qpcR)
library(ggplot2)

data_cells <- readRDS("SC_topo_NaCl_anno.rds")
dim(data_cells)
# 16566 35932

# fix mouse ID
data_cells$mouse_label <- data_cells$SampleID
data_cells$mouse_label <- gsub(x = data_cells$mouse_label, pattern = "Mouse25T", replacement = "Mouse4_HSD", fixed = TRUE)

data_cells$mouse_label <- gsub(x = data_cells$mouse_label, pattern = "Mouse32T", replacement = "Mouse3_HSD", fixed = TRUE)
  
data_cells$mouse_label <- gsub(x = data_cells$mouse_label, pattern = "Mouse4NT", replacement = "Mouse1_NSD", fixed = TRUE)
  
data_cells$mouse_label <- gsub(x = data_cells$mouse_label, pattern = "Mouse5NT", replacement = "Mouse2_NSD", fixed = TRUE)  

# frequencies of each cluster in each mouse
library(RColorBrewer)
library(ggplot2)

tab <- prop.table(table(data_cells$desc_0.5,data_cells$mouse_label), margin = 2)

#write.table(tab, "raw_data_FIG2_D.txt", sep="\t", col.name = TRUE, quote = FALSE)
#tab <- read.table(file="raw_data_FIG2_D.txt", header= TRUE, stringsAsFactors=FALSE, sep="\t")

tab <- t(tab)

n <- nrow(tab)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

pdf("FIG2_D.pdf", width = 10, height = 10)
barplot(tab, beside=TRUE, ylim=c(0, max(tab) + 0.1), col= sample(col_vector, n), legend.text=TRUE, bty = "n")
dev.off()

data_cells <- NormalizeData(data_cells, normalization.method = "LogNormalize", scale.factor = 10000)

data_cells <- FindVariableFeatures(data_cells, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(data_cells), 10)

all.genes <- rownames(data_cells)
data_cells <- ScaleData(data_cells, features = all.genes)

Idents(object = data_cells) <- "desc_0.5"
data_cells$desc05_treat <- paste(Idents(data_cells), data_cells$Treatment, sep = "_")
Idents(object = data_cells) <- "desc05_treat"

# load original umap coordinates
coord_original <- read.table(file="raw_data_FIG2_C.txt", header= TRUE, stringsAsFactors=FALSE, sep="\t", row.names = 1)

data_cells <- AddMetaData(data_cells, metadata = coord_original)
head(data_cells@meta.data, 5)

# make a new umap to have the slot

data_cells <- RunPCA(data_cells, features = VariableFeatures(object = data_cells))

data_cells <- RunUMAP(data_cells, dims = 1:20)

# overwrite the umap slot with the original coord
coord <- data_cells@meta.data[,35:36]
colnames(coord) <- c("umaporiginal_1","umaporiginal_2")
head(coord)

dim(coord)
# 35932     2

coord <- as.matrix(coord)	
data_cells[["umap"]] <- CreateDimReducObject(embeddings = coord, key = "umap_original_", assay = DefaultAssay(data_cells))

head(data_cells[["umap"]]@cell.embeddings)

Idents(data_cells) <- "desc_0.5"

# UMAP split by treatment
pdf("FIG2_C.pdf", width = 25)
DimPlot(data_cells, reduction = "umap", split.by="Treatment", label = TRUE)
dev.off()


##########################################
####### marker cluster 4 HSD vs NSD

C4_cells_markers <- FindMarkers(data_cells, ident.1 = "4_Treated", ident.2 = "4_NotTreated", verbose = FALSE, logfc.threshold= 0)

C4_cells_markers_1 <- FindMarkers(data_cells, ident.1 = "4_Treated", ident.2 = "4_NotTreated", verbose = FALSE)

C4_cells_markers_1_fdr001 <- C4_cells_markers_1[C4_cells_markers_1$p_val_adj <= 0.01,]
write.table(C4_cells_markers_1_fdr001, "DEGs_Cluster4_treatVSNOTtreat_markers_fdr001.txt", sep="\t")


##########################################
##########################################

#VOLCANO plot
# FDR 0.05

volcano_data <- C4_cells_markers

#write.table(volcano_data, "raw_data_FIG2_E.txt", sep="\t", col.name = TRUE, quote = FALSE)
#tab_2E <- read.table(file="raw_data_FIG2_E.txt", header= TRUE, stringsAsFactors=FALSE, sep="\t")

not_sign <- C4_cells_markers$p_val_adj > 0.05
up <- C4_cells_markers$avg_log2FC >= 0.5 & C4_cells_markers$p_val_adj <= 0.05
down <- C4_cells_markers$avg_log2FC <= -0.5 & C4_cells_markers$p_val_adj <= 0.05

gene_selected_c4 <- c("Irf8", "Prf1", "Stat3", "Id2", "Nfkb1", "Csf1", "Klrk1", "Klrd1", "Oas1a", "Eomes", "Prdm1", "Maf", "Cxcr3", "Jun", "Klf2", "Gzmk", "Gzmf", "Gzme", "Gzmd", "Gzmc", "Ifitm1", "Tnfrsf9")
length(gene_selected_c4)
# 22

results <- volcano_data
results$gene <- rownames(results)
results_plot <- results[gene_selected_c4,]
dim(results_plot)
# 22  6

pdf("FIG2_E.pdf")
plot(volcano_data$avg_log2FC, -1*log10(volcano_data$p_val_adj), main="Volcano plot Cluster 4", xlab="avg log2FC", ylab="-log10(Adjusted P-Value)", xlim = c(-(round(max(abs(volcano_data$avg_log2FC)))),round(max(abs(volcano_data$avg_log2FC)))))
points(volcano_data$avg_log2FC[not_sign], -1*log10(volcano_data$p_val_adj[not_sign]), col="gray87")
points(volcano_data$avg_log2FC[up], -1*log10(volcano_data$p_val_adj[up]), col="red", pch = 20)
points(volcano_data$avg_log2FC[down], -1*log10(volcano_data$p_val_adj[down]), col="blue", pch = 20)
abline(h=c(1.30103), col="black")
abline(v=c(-0.5, 0.5), col="black")
text(results_plot$avg_log2FC, -1*log10(results_plot$p_val_adj),labels = results_plot$gene, col="black",cex=0.5,pos=3)
dev.off()



##########################################
####### rnk GSEA


C4_cells_markers_rnk <- FindMarkers(object = data_cells, ident.1= "4_Treated", ident.2= "4_NotTreated", logfc.threshold = 0, min.pct = 0)

C4_cells_rank <- C4_cells_markers_rnk[order(as.numeric(C4_cells_markers_rnk$avg_log2FC),decreasing = TRUE),]


C4_cells_rank <-cbind(rownames(C4_cells_rank),C4_cells_rank$avg_log2FC)
colnames(C4_cells_rank) <- c("GeneName","rank")

C4_cells_rank <- as.data.frame(C4_cells_rank)

write.table(C4_cells_rank, file="raw_data_FIG2_F.rnk", sep="\t", row.names=FALSE, col.name = TRUE, quote = FALSE)


##########################################
### subcluster C4
# filter only CD8 (cluster 4)
head(data_cells@meta.data, 5)

CD8_cells <- subset(data_cells, subset = anno_desc_0.5 == "CD8_Tcells")
dim(CD8_cells)
# 16566  2508

all.genes <- rownames(CD8_cells)
CD8_cells <- FindVariableFeatures(CD8_cells, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(CD8_cells), 10)

CD8_cells <- ScaleData(CD8_cells, features = all.genes)

CD8_cells <- RunPCA(CD8_cells, features = VariableFeatures(object = CD8_cells))

CD8_cells <- JackStraw(CD8_cells, num.replicate = 100)
CD8_cells <- ScoreJackStraw(CD8_cells, dims = 1:20)

CD8_cells <- FindNeighbors(CD8_cells, reduction = "pca", dims = 1:20)

##### find different resolutions

CD8_cells <- FindClusters(CD8_cells, resolution = seq(from=0, to=2, by=0.1), print.output = 0, save.SNN = T)

head(CD8_cells@meta.data, 5)

library(clustree)

Idents(CD8_cells) <- "RNA_snn_res.0.1"

#########################
# remove cluster 4 della res 0.1 (gamma/delta t cells) 
#########################

table(CD8_cells$RNA_snn_res.0.1)
#   0   1   2   3   4 
# 810 656 589 299 154 

tcells_clean <- subset(CD8_cells, subset = RNA_snn_res.0.1 == 4, invert = TRUE)
table(tcells_clean$RNA_snn_res.0.1)
dim(tcells_clean)
# 16566  2354

all.genes <- rownames(tcells_clean)
tcells_clean <- FindVariableFeatures(tcells_clean, selection.method = "vst", nfeatures = 2000)
tcells_clean <- ScaleData(tcells_clean, features = all.genes)

head(tcells_clean@meta.data, 5)

Idents(tcells_clean) <- "RNA_snn_res.0.1"

##### UMAP

tcells_clean <- RunUMAP(tcells_clean, dims = 1:20)

tab_2H <- cbind(tcells_clean@reductions$umap@cell.embeddings, tcells_clean$RNA_snn_res.0.1)
colnames(tab_2H) <- c("dim1","dim2","cluster")

#write.table(tab_2H, "raw_data_FIG2_H.txt", sep="\t", col.name = TRUE, quote = FALSE)

pdf("FIG2_H.pdf")
DimPlot(tcells_clean, reduction = "umap", group.by="RNA_snn_res.0.1", label = TRUE)
dev.off()

# frequencies of each cluster in each mouse
library(RColorBrewer)
library(ggplot2)
tab <- prop.table(table(tcells_clean$RNA_snn_res.0.1,tcells_clean$SampleID), margin = 2)
n <- nrow(tab)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

raw_2I <- tab
colnames(raw_2I) <- c("Mouse4_HSD","Mouse3_HSD","Mouse1_NSD","Mouse2_NSD")
raw_2I <- raw_2I[-5,]
raw_2I <- raw_2I[,c("Mouse1_NSD","Mouse2_NSD","Mouse4_HSD","Mouse3_HSD")]

#write.table(raw_2I, "raw_data_FIG2_I.txt", sep="\t", col.name = TRUE, quote = FALSE)
#tab <- read.table(file="raw_data_FIG2_I.txt", header= TRUE, stringsAsFactors=FALSE, sep="\t")
#tab <- as.matrix(tab)

pdf("FIG2_I.pdf", width = 10, height = 10)
barplot(tab, xlim=c(0, ncol(tab) + 3), col= sample(col_vector, n), legend.text=TRUE, args.legend=list(x=ncol(tab) + 3, y=max(colSums(tab)), bty = "n"))
dev.off()


markers <- c("Tox","Tigit","Ctla4","Pdcd1","Havcr2","Entpd1","Mki67","Gzmb","Gzma","Gzmk","Zfp683","Cx3cr1","Tbx21","Klrg1","Ifng","Klrb1c","Nkg7","Runx3","Cd27","Il7r","Ccr7","Sell","Tcf7","Cd69","Itgae")


pdf("FIG2_J.pdf", width=13)
DotPlot(tcells_clean, features = markers) + RotatedAxis()
dev.off()

gg_data <- DotPlot(tcells_clean, features = markers)
#write.table(gg_data$data, "raw_data_FIG2_J.txt", sep="\t", col.name = TRUE, quote = FALSE)

tcells_clean_markers <- FindAllMarkers(tcells_clean, only.pos = TRUE, return.thresh = 0.05)
write.table(tcells_clean_markers, "CD8_cells_clean_res01_markers.txt", sep="\t")

top20 <- tcells_clean_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

tcells_clean_markers_fdr001 <- tcells_clean_markers[tcells_clean_markers$p_val_adj <= 0.01,]
write.table(tcells_clean_markers_fdr001, "CD8_cells_clean_res01_markers_fdr001.txt", sep="\t")
