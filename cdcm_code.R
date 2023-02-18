## Download and load packages needed for this analysis

library(Matrix)
library(Seurat)
library(DropletUtils)
library(ggplot2)
library(ggthemes)
library(scales)
library(dplyr)
library(deMULTIplex)
library(patchwork)
library(clipr)
library(reshape2)
library(stringr)
library(monocle)
library(tradeSeq)
library(pheatmap)
library(biomaRt)
library(sctransform)
library(reshape2)
library(grid)
library(monocle)
library(reshape2)

#For each dataset, load in, and append a label to the column name so you know which group it is from
data1 = readMM("~/Desktop/honglab/seqexpt/filtered_feature_bc_matrix_CC2/matrix.mtx.cc2.gz")
rownames(data1) = read.table("~/Desktop/honglab/seqexpt/filtered_feature_bc_matrix_CC2/features.tsv.cc2.gz", as.is = TRUE)$V1
colnames(data1) = read.table("~/Desktop/honglab/seqexpt/filtered_feature_bc_matrix_CC2/barcodes.tsv.cc2.gz", as.is = TRUE)$V1
colnames(data1) = paste(colnames(data1), "control", sep = "_")


data2 = readMM("~/Desktop/honglab/seqexpt/filtered_feature_bc_matrix_LVIP/matrix.mtx.lvip.gz")
rownames(data2) = read.table("~/Desktop/honglab/seqexpt/filtered_feature_bc_matrix_LVIP/features.tsv.lvip.gz", as.is = TRUE)$V1
colnames(data2) = read.table("~/Desktop/honglab/seqexpt/filtered_feature_bc_matrix_LVIP/barcodes.tsv.lvip.gz", as.is = TRUE)$V1
colnames(data2) = paste(colnames(data2), "cDCM", sep = "_")

#combine the two datasets, then remove the intermediates
data = cbind(data1, data2)
rm(data1, data2)

#make the phenotype table
phenotype = data.frame(row.names = colnames(data), group = sapply(strsplit(colnames(data), "_"), "[[", 2))

#Rename genes using biomaRt
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list = getBM(filters = "ensembl_gene_id", attributes = c("hgnc_symbol", "ensembl_gene_id"), values = rownames(data), mart = mart, useCache = FALSE) #have to put useCache = FALSE because my bioconductor and R are behind; probably going to need to update soon
G_list = G_list[!(duplicated(G_list$hgnc_symbol)) & !(is.na(G_list$hgnc_symbol)) & !(G_list$hgnc_symbol == ""), ]
data = as.data.frame(as.matrix(data)) #this step requires a lot of ram
data$ensembl_gene_id = rownames(data)
data_named = merge(data, G_list, by = "ensembl_gene_id", all.x = FALSE)
rownames(data_named) = data_named$hgnc_symbol
data_named$hgnc_symbol = NULL
data_named$ensembl_gene_id = NULL
data = data_named
rm(data_named, G_list, mart)
data = Matrix(as.matrix(data), sparse = TRUE)


### Begin seurat vignette ( we used version 3.2.3)
pbmc <- CreateSeuratObject(counts = data, project = "DCM", min.cells = 3, min.features = 200, meta.data = phenotype[1:21868, , drop = FALSE])
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
### We customized QC metrics as cardiomyocytes in 10x chips have higher mitochondria counts than other cell types
pbmc <- subset(pbmc, subset = nFeature_RNA > 2500 & nFeature_RNA < 8750 & percent.mt < 35)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.2)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)
table(pbmc$group)


#### filter non-cardiac and cell cycle clusters
pbmc.clean = pbmc[, !pbmc$seurat_clusters %in% c(1, 3, 5)] 
pbmc.clean = RunPCA(pbmc.clean, verbose = FALSE)
pbmc.clean = RunUMAP(pbmc.clean, dims = 1:30)
pbmc.clean = FindNeighbors(pbmc.clean, dims = 1:15)
pbmc.clean = FindClusters(pbmc.clean, resolution = .15) 
DimPlot(pbmc.clean, label = TRUE, pt.size = 1)
DimPlot(pbmc.clean, reduction = "umap", label = TRUE)
DimPlot(pbmc.clean, group.by = "group")

### load in present cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

### visualize clusters undergoing cell cycle
marrow <- CellCycleScoring(pbmc.clean, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(marrow)
DimPlot(marrow, group.by = "seurat_clusters")

### remove cycling cluster
pbmc.clean.1 = pbmc.clean[, !pbmc.clean$seurat_clusters %in% c(3)] 
pbmc.clean.1 = RunPCA(pbmc.clean.1, verbose = FALSE)
pbmc.clean.1 = RunUMAP(pbmc.clean.1, dims = 1:30)
pbmc.clean.1 = FindNeighbors(pbmc.clean.1, dims = 1:15)
pbmc.clean.1 = FindClusters(pbmc.clean.1, resolution = .15) 
DimPlot(pbmc.clean.1, label = TRUE, pt.size = 1)
DimPlot(pbmc.clean.1, reduction = "umap", label = TRUE)
DimPlot(pbmc.clean.1, group.by = "group")

#### integrate datasets
options(future.globals.maxSize = 4000 * 1024^2)
pbmc.list = SplitObject(pbmc.clean.1, split.by = "group")
for (i in 1:length(pbmc.list)) {
  pbmc.list[[i]] = NormalizeData(pbmc.list[[i]], verbose = FALSE)
  pbmc.list[[i]] = FindVariableFeatures(pbmc.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}
for (i in 1:length(pbmc.list)) {
  pbmc.list[[i]] = SCTransform(pbmc.list[[i]], verbose = FALSE)
}
pbmc.features = SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 3000)
pbmc.list = PrepSCTIntegration(object.list = pbmc.list, anchor.features = pbmc.features, 
                               verbose = FALSE)
pbmc.anchors = FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT", 
                                      anchor.features = pbmc.features, verbose = FALSE)
pbmc.integrated = IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT", 
                                verbose = FALSE)
pbmc.integrated = RunPCA(pbmc.integrated, verbose = FALSE)
pbmc.integrated = RunUMAP(pbmc.integrated, dims = 1:30)
DimPlot(pbmc.integrated, group.by = "group")
pbmc.integrated = FindNeighbors(pbmc.integrated, dims = 1:30)
pbmc.integrated = FindClusters(pbmc.integrated, resolution = .20)
DimPlot(pbmc.integrated, label = TRUE)
rm(pbmc.anchors) #clear some memory
DimPlot(pbmc.integrated, group.by = "group")
DimPlot(pbmc.integrated, label = "TRUE")

### filter cluster with artefactual apoptosis/cell death signatures ###
pbmc.integrated.clean = pbmc.integrated[, !pbmc.integrated$seurat_clusters %in% c(5)] 
pbmc.integrated.clean = RunPCA(pbmc.integrated.clean, verbose = FALSE)
pbmc.integrated.clean = RunUMAP(pbmc.integrated.clean, dims = 1:30)
all.genes = rownames(pbmc.integrated.clean)
pbmc.integrated.clean <- ScaleData(pbmc.integrated.clean, features = all.genes)
pbmc.integrated.clean = FindNeighbors(pbmc.integrated.clean, dims = 1:15)
pbmc.integrated.clean = FindClusters(pbmc.integrated.clean, resolution = .15) 
DimPlot(pbmc.integrated.clean, label = TRUE, pt.size = 1)
DimPlot(pbmc.integrated.clean, group.by = "group")
DimPlot(pbmc.integrated.clean, label = "TRUE")

DefaultAssay(pbmc.integrated.clean) = "RNA"

### now that we have filtered appropriately we can identify differentially expressed genes between the two groups
lvip.markers <- FindMarkers(pbmc.integrated.clean, group.by = "group", ident.1 = "cDCM", logfc.threshold = 0.25, only.pos = TRUE)
cc2.markers <- FindMarkers(pbmc.integrated.clean, group.by = "group", ident.1 = "control", logfc.threshold = 0.25, only.pos = TRUE)

### generating entropy score. for full explanation on calculating entropy score please see cited paper.

### subset only CMS for entropy 
goodcms_control_entropy = data1[, colnames(data1) %in% colnames(pbmc.integrated.clean)]
goodcms_cDCM_entropy = data2[, colnames(data2) %in% colnames(pbmc.integrated.clean)]

gene_control = rename_genes(goodcms_control_entropy, species = "human")
gene_cDCM = rename_genes(goodcms_cDCM_entropy, species = "human")

entropy_gene_control = master_entropy(gene_control, species = "human")
entropy_gene_cDCM = master_entropy(gene_cDCM, species = "human")

a = as.data.frame(entropy_gene_control)
b = as.data.frame(entropy_gene_cDCM)

colnames(a) = "score"
colnames(b) = "score"
vec = "control"
vec1 = "lvip"
a1 = a 
a1$condition = vec
b1 = b
b1$condition = vec1
c = rbind(a1,b1)
ggplot(c, aes(x= condition, y = score)) + geom_boxplot()

d <-ggplot(c, aes(x=condition, y=score, fill=condition)) + geom_boxplot()

d + scale_fill_manual(values = c("thistle1", "lightskyblue1")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
 d

#### Generating heatmap based on candidate genes
DoHeatmap(pbmc.integrated.clean, group.by = "group", features = c("MYH6", "TNNI3", "NKX2-5", "TNNT2", "TTN", "MYH7", "MYL2", "MYL7", "TNNI1",  "MYOM1", "MT-CO1", "MT-CO2", "MT-CO3", "MT-ND2", "MT-ND5", "MT-ND6", "MT-CYB", "MT-ND4", "COX6C", "COX6A2"), group.colors = c("thistle1", "lightskyblue1"), raster = FALSE) + scale_fill_gradientn(colors = c("royalblue", "white", "indianred1"))


### after inputting DE genes on geneontology.org, we use this to make our dot plot with the generated GO terms
data1 <- read.csv("~/Desktop/honglab/go_actn2_iso2.csv", header=TRUE, stringsAsFactors = FALSE)
library("ggplot2")
S1<- ggplot(data1, aes(x=DE, y=go_term, size=fe, color=pval)) +
  geom_point(alpha = .8)  +
  # facet_grid(cmko ~ .) +
  theme_classic()
S1 = S1 + theme(panel.grid.major.x=element_line(size=.5,colour="grey88"), panel.grid.major.y=element_line(size=0.5,colour="grey88"))
S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab")
ggsave("goterms_down.svg",plot=S1, height = 6, width = 7)
