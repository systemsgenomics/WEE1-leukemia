# mouse PDX data, checked  here human cells
#cell ranger GRCh38-2020-A #Pipeline Version	cellranger-6.0.1

# MLL7 KMT2A:AFF1 PDX cells treated with wee1i 24h + 4h

### Seurat package analysis
### http://satijalab.org/seurat/

library(Seurat)
library(dplyr)
library(Matrix)
#library(gdata)


#MLL7 Treated -----
matrix_dir = "p2_12/outs/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")

mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,header = FALSE,stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,header = FALSE,stringsAsFactors = FALSE)


head(barcode.names)

head(feature.names)

colnames(mat) = barcode.names$V1
rownames(mat)= feature.names[,2]

head(mat)


MLL7tcounts=colSums(mat)
summary(MLL7tcounts)


MLL7t <- CreateSeuratObject(counts = mat, project = "PDX", min.cells = 2, min.features = 10)
MLL7t=AddMetaData(MLL7t, rep("MLL7treated",nrow(MLL7t[[]])),"trt")
head(MLL7t[[]])

rm(mat)
### ----------------------------------------

##NO  MERGING
MLL7tfilt=MLL7t

MLL7tfilt[["percent.mt"]] <- PercentageFeatureSet(MLL7tfilt, pattern = "^MT-")

summary(MLL7tfilt$percent.mt)

pdf("VlnPlot_nFeature_MLL7treated.pdf")
VlnPlot(MLL7tfilt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4)
VlnPlot(MLL7tfilt, features = "percent.mt")
VlnPlot(MLL7tfilt, features = "nFeature_RNA")
dev.off()
#

MLL7tfilt <- subset(MLL7tfilt, subset = nFeature_RNA > 1500 & nFeature_RNA < 6000 & percent.mt < 10)

table(MLL7tfilt$trt)

## NORMALIZE. Values stored in seurat[["RNA"]]@data
MLL7tfilt <- NormalizeData(MLL7tfilt, normalization.method = "LogNormalize", scale.factor = 10000)

## Identification of highly variable features (feature selection)
MLL7tfilt <- FindVariableFeatures(MLL7tfilt, selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(MLL7tfilt))
top30 <- head(VariableFeatures(MLL7tfilt), 30)
# plot variable features with and without labels
pdf("top30varGenes_MLL7tfilt.pdf")
plot1=VariableFeaturePlot(MLL7tfilt)
plot2=LabelPoints(plot = plot1, points = top30, repel = TRUE)
plot2
dev.off()




## PLOT CELL CYCLE GENES
## list of cell cycle markers, from Tirosh et al, 2015
## S -> G2 -> M -> G1 -> <- G0 --> S
## Cell cycle arrest = G1 -> G0
## M = Mitosis (nuclear division) + cytokinesis (division); S = synthesis (DNA replication)
s.genes <- c("MCM5","PCNA","TYMS", "FEN1","MCM2" ,"MCM4","RRM1","UNG", "GINS2","MCM6" ,"CDCA7" ,"DTL",
               "PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1", "GMNN" ,"WDR76" ,"SLBP","CCNE2",
               "UBR7","POLD3","MSH2", "ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM",
               "CASP8AP2","USP1","CLSPN","POLA1" , "CHAF1B","BRIP1","E2F8")
g2m.genes <- c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5",
                 "TPX2","TOP2A","NDC80", "CKS2", "NUF2","CKS1B","MKI67","TMPO","CENPF" ,"TACC3","FAM64A" ,"SMC4",
                 "CCNB2", "CKAP2L","CKAP2", "AURKB" ,"BUB1" ,"KIF11","ANP32E","TUBB4B" , "GTSE1" ,"KIF20B","HJURP","CDCA3",
                 "HN1" ,"CDC20","TTK" , "CDC25C","KIF2C","RANGAP1", "NCAPD2" ,"DLGAP5","CDCA2" , "CDCA8" ,"ECT2" ,"KIF23",
                 "HMMR", "AURKA","PSRC1","ANLN","LBR", "CKAP5" , "CENPE" , "CTCF" ,"NEK2" ,"G2E3" , "GAS2L3","CBX5","CENPA")

MLL7tfilt <- CellCycleScoring(object = MLL7tfilt, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)

pdf("cellcycleGenes_MLL7tfilt.pdf")
FeatureScatter(MLL7tfilt, "S.Score", "nCount_RNA")
dev.off()


## Running a PCA on cell cycle genes
all.genes <- rownames(MLL7tfilt)
MLL7tfilt <- ScaleData(MLL7tfilt, features = all.genes)
MLL7tfilt <- RunPCA(MLL7tfilt, features = c(s.genes, g2m.genes))
pdf("PCA_cellCycleGenes_MLL7tfilt.pdf")
DimPlot(MLL7tfilt,reduction = "pca")
DimPlot(MLL7tfilt,reduction = "pca", group.by = 'trt')
dev.off()


table(paste0(MLL7tfilt$trt,MLL7tfilt$Phase))

MLL7tfilt=AddMetaData(MLL7tfilt, paste(MLL7tfilt$trt,MLL7tfilt$Phase, sep=";"),"group")

## Perform linear dimensional reduction
pdf("PCA_varGenes_MLL7tfilt.pdf")
MLL7tfilt <- RunPCA(MLL7tfilt, features = VariableFeatures(object = MLL7tfilt), nfeatures.print = 5)
VizDimLoadings(MLL7tfilt, dims = 1:2, reduction = "pca")
DimPlot(MLL7tfilt, reduction = "pca", group.by = 'group', dims = c(1, 2),pt.size = 0.3)
DimPlot(MLL7tfilt, reduction = "pca", group.by = 'Phase')
DimHeatmap(MLL7tfilt, dims = 1:6, cells = 500, balanced = TRUE)

## Determine the ‘dimensionality’ of the dataset - overcome technical noise
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
MLL7tfilt <- JackStraw(MLL7tfilt, num.replicate = 50, dims=30)
MLL7tfilt <- ScoreJackStraw(MLL7tfilt, dims = 1:30)
JackStrawPlot(MLL7tfilt, dims = 1:30)
ElbowPlot(MLL7tfilt)
dev.off()

## Cluster the cells
MLL7tfilt <- FindNeighbors(MLL7tfilt, dims = 1:10)
# was 0.1 before, lets increase a bit
MLL7tfilt <- FindClusters(MLL7tfilt, resolution = 0.5)  # resolution - reducing the number of clusters found in UMAP


## Run non-linear dimensional reduction (UMAP/tSNE)
MLL7tfilt <- RunUMAP(MLL7tfilt, reduction = "pca", dims = 1:6, n_neighbors = 10)
pdf("UMAPvarGenes_MLL7tfilt.pdf")#has 0.5 resolution, best now
DimPlot(MLL7tfilt, reduction = "umap", label=T)
DimPlot(MLL7tfilt, reduction = "umap", group.by = 'group', pt.size = 0.3)
DimPlot(MLL7tfilt, reduction = "umap", group.by = 'trt', pt.size = 0.3)
## Highlight by Cell Cycle Phase
DimPlot(MLL7tfilt, reduction = "umap", group.by = 'Phase')
dev.off()


#save the Seurat object here
saveRDS(MLL7tfilt, file = "MLL7treated.filt.rds")

##check cell amounts
table(paste0(MLL7tfilt$trt))


table(paste0(MLL7tfilt$trt,MLL7tfilt$Phase))

table(paste0(MLL7tfilt$seurat_clusters))

pdf("QC_UMAPvarGenes_MLL7treated.filt.pdf")
FeaturePlot(MLL7tfilt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()

##Make Findmarkers per cluster from single samples

markers <- FindAllMarkers(MLL7tfilt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = 'RNA')
#save(markers, file="markers.RData")
#need dplyr package to install

markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top1000 <- markers %>% group_by(cluster) %>% top_n(n = 1000, wt = avg_logFC)
class(top1000)
write.table(top1000, file="top1000clusterMarkers_MLL7treated.txt", sep="\t", row.names=F, quote=F)


## Examine genes of interest

#pdf("BCL6_MLL7tfilt.pdf")
#FeaturePlot(MLL7tfilt, features = "BCL6")
#dev.off()



#########################
