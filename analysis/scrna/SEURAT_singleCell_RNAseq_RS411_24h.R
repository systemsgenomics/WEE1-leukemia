## COMBINE DMSO & WEE1i SAMPLES (24 h)

### Seurat package analysis

library(Seurat)
library(dplyr)
library(Matrix)


## sample 8, DMSO ----------------------
matrix_dir = "cellranger/S8/outs/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")

mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,header = FALSE,stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,header = FALSE,stringsAsFactors = FALSE)

cell_meta=read.table("cellranger/S8/outs/analysis/gem_classification.csv", sep=",", header=T, stringsAsFactors=F)

## Data has both human and mouse gene counts, and human and mouse cells
hg19c=cell_meta$call=="hg19"

hg19g=1:nrow(feature.names)%in%grep("hg19", feature.names[,1])

hg19counts=colSums(mat[hg19g,hg19c])
summary(hg19counts)


# sanity check, calc mouse gene counts from human cells
mm10g=1:nrow(feature.names)%in%grep("mm10", feature.names[,1])
mm10counts=colSums(mat[mm10g,hg19c])
summary(mm10counts)


feature.names[,2]=sub("hg19_","",feature.names[,2])
feature.names[,1]=sub("hg19_","",feature.names[,1])

## matching to symbols not unique, fix here
feature.names[duplicated(feature.names[,2]),2]=paste(feature.names[duplicated(feature.names[,2]),2],feature.names[duplicated(feature.names[,2]),1], sep=":")

colnames(mat) = barcode.names$V1
rownames(mat)= feature.names[,2]

mat=mat[hg19g,hg19c]


D <- CreateSeuratObject(counts = mat, project = "RS411",  min.features = 10)
D=AddMetaData(D, rep("DMSO",nrow(D[[]])),"trt")
D=AddMetaData(D, mm10counts,"mm10counts")
head(D[[]])
### ----------------------------------------

## sample 6, WEE1i ----------------------
matrix_dir = "cellranger/S6/outs/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")

mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,header = FALSE,stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,header = FALSE,stringsAsFactors = FALSE)

cell_meta=read.table("//cellranger/S6/outs/analysis/gem_classification.csv", sep=",", header=T, stringsAsFactors=F)

## Data has both human and mouse gene counts, and humand and mouse cells
hg19c=cell_meta$call=="hg19"

hg19g=1:nrow(feature.names)%in%grep("hg19", feature.names[,1])

hg19counts=colSums(mat[hg19g,hg19c])
summary(hg19counts)


# sanity check, calc mouse gene counts from human cells
mm10g=1:nrow(feature.names)%in%grep("mm10", feature.names[,1])
mm10counts=colSums(mat[mm10g,hg19c])
summary(mm10counts)

feature.names[,2]=sub("hg19_","",feature.names[,2])
feature.names[,1]=sub("hg19_","",feature.names[,1])

## matching to symbols not unique, fix here
feature.names[duplicated(feature.names[,2]),2]=paste(feature.names[duplicated(feature.names[,2]),2],feature.names[duplicated(feature.names[,2]),1], sep=":")

colnames(mat) = barcode.names$V1
rownames(mat)= feature.names[,2]

mat=mat[hg19g,hg19c]


W <- CreateSeuratObject(counts = mat, project = "RS411",  min.features = 10)
W=AddMetaData(W, rep("WEE1i",nrow(W[[]])),"trt")
W=AddMetaData(W, mm10counts,"mm10counts")
head(W[[]])
rm(mat)
### ----------------------------------------

## MERGE
RS=merge(D,W, min.cells=10)

## stash QC stats
RS[["percent.mt"]] <- PercentageFeatureSet(RS, pattern = "^MT-")

pdf("VlnPlot_nFeature.pdf")
VlnPlot(RS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","mm10counts"), ncol = 4)
VlnPlot(RS, features = "percent.mt", split.by = "trt")
VlnPlot(RS, features = "mm10counts", split.by = "trt")
dev.off()
#

## FILTER by QC
RS.filt <- subset(RS, subset = nFeature_RNA > 2000 & nFeature_RNA < 6000 & percent.mt < 20 & mm10counts < 500)
table(RS.filt$trt)

## NORMALIZE. Values stored in seurat[["RNA"]]@data
RS.filt <- NormalizeData(RS.filt, normalization.method = "LogNormalize", scale.factor = 10000)

## Identification of highly variable features (feature selection)
RS.filt <- FindVariableFeatures(RS.filt, selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(RS.filt))
top30 <- head(VariableFeatures(RS.filt), 30)
# plot variable features with and without labels
pdf("top30varGenes.pdf")
plot1=VariableFeaturePlot(RS.filt)
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

RS.filt <- CellCycleScoring(object = RS.filt, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)

pdf("cellcycleGenes.pdf")
FeatureScatter(RS.filt, "S.Score", "nCount_RNA")
dev.off()


## Running a PCA on cell cycle genes
all.genes <- rownames(RS.filt)
RS.filt <- ScaleData(RS.filt, features = all.genes)
RS.filt <- RunPCA(RS.filt, features = c(s.genes, g2m.genes))
pdf("PCA_cellCycleGenes.pdf")
DimPlot(RS.filt,reduction = "pca")
DimPlot(RS.filt,reduction = "pca", group.by = 'trt')
dev.off()

table(paste0(RS.filt$trt,RS.filt$Phase))

RS.filt=AddMetaData(RS.filt, paste(RS.filt$trt,RS.filt$Phase, sep=";"),"group")

## Perform linear dimensional reduction
pdf("PCA_varGenes.pdf")
RS.filt <- RunPCA(RS.filt, features = VariableFeatures(object = RS.filt), nfeatures.print = 5)
VizDimLoadings(RS.filt, dims = 1:2, reduction = "pca")
DimPlot(RS.filt, reduction = "pca", group.by = 'group', dims = c(1, 2),pt.size = 0.3)


DimPlot(RS.filt, reduction = "pca", group.by = 'Phase')
DimHeatmap(RS.filt, dims = 1:6, cells = 500, balanced = TRUE)

## Determine the ‘dimensionality’ of the dataset - overcome technical noise
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
RS.filt <- JackStraw(RS.filt, num.replicate = 50, dims=30)
RS.filt <- ScoreJackStraw(RS.filt, dims = 1:30)
JackStrawPlot(RS.filt, dims = 1:30)
## heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one
ElbowPlot(RS.filt)
dev.off()

## Cluster the cells
RS.filt <- FindNeighbors(RS.filt, dims = 1:10)
RS.filt <- FindClusters(RS.filt, resolution = 0.5)  # resolution - reducing the number of clusters found in UMAP

## Run non-linear dimensional reduction (UMAP/tSNE)
RS.filt <- RunUMAP(RS.filt, reduction = "pca", dims = 1:6, n_neighbors = 10)

pdf("UMAPvarGenes.pdf")
DimPlot(RS.filt, reduction = "umap", label=T)
DimPlot(RS.filt, reduction = "umap", group.by = 'group', pt.size = 0.3)
## Highlight by Cell Cycle Phase
DimPlot(RS.filt, reduction = "umap", group.by = 'Phase')
dev.off()

saveRDS(RS.filt, file = "RSfilt.rds")

pdf("QC_UMAPvarGenes_old.pdf")
FeaturePlot(RS.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","mm10counts"))
dev.off()

## Examine genes of interest
#readRDS("RSfilt.rds")
## Featureplots
#FeaturePlot(RS.filt, features = c("TOP2A", "MDM2"))
#dev.off()


## Finding differentially expressed features (cluster biomarkers)

# find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(RS.filt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = 'RNA')
markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#pdf("markerGenes.pdf")
#DoHeatmap(RS.filt, features = top10$gene) + NoLegend()
#dev.off()

