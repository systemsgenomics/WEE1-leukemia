## Nalm6

library(Seurat)
library(dplyr)
library(Matrix)
#library(gdata)

## sample 2, DMSO ----------------------
matrix_dir = "cellranger/S2/outs/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")

mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,header = FALSE,stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,header = FALSE,stringsAsFactors = FALSE)

cell_meta=read.table("cellranger/S2/outs/analysis/gem_classification.csv", sep=",", header=T, stringsAsFactors=F)


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


D <- CreateSeuratObject(counts = mat, project = "RS411",  min.features = 10)
D=AddMetaData(D, rep("DMSO",nrow(D[[]])),"trt")
D=AddMetaData(D, mm10counts,"mm10counts")
head(D[[]])
### ----------------------------------------

## sample 3, WEE1i ----------------------
matrix_dir = "cellranger/S3/outs/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")

mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,header = FALSE,stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,header = FALSE,stringsAsFactors = FALSE)

cell_meta=read.table("cellranger/S3/outs/analysis/gem_classification.csv", sep=",", header=T, stringsAsFactors=F)

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


## MERGE: https://satijalab.org/seurat/v3.1/merge_vignette.html
## lets first just use the regular nalm6
N=merge(D,W, min.cells=10)


## The [[ operator can add columns to object metadata. This is a great place to stash QC stats
N[["percent.mt"]] <- PercentageFeatureSet(N, pattern = "^MT-")

pdf("VlnPlot_nFeature.pdf")
VlnPlot(N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","mm10counts"), ncol = 4)
VlnPlot(N, features = "percent.mt", split.by = "trt")
VlnPlot(N, features = "mm10counts", split.by = "trt")
dev.off()
#

## SUBSET = FILTER & check that better
N.filt <- subset(N, subset = nFeature_RNA > 2000 & nFeature_RNA < 6000 & percent.mt < 20 & mm10counts < 500)

table(N.filt$trt)


## NORMALIZE. Values stored in seurat[["RNA"]]@data
N.filt <- NormalizeData(N.filt, normalization.method = "LogNormalize", scale.factor = 10000)
## Identification of highly variable features (feature selection)
N.filt <- FindVariableFeatures(N.filt, selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(N.filt))
top30 <- head(VariableFeatures(N.filt), 30)
# plot variable features with and without labels
pdf("top30varGenes.pdf")
plot1=VariableFeaturePlot(N.filt)
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

N.filt <- CellCycleScoring(object = N.filt, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)

pdf("cellcycleGenes.pdf")
FeatureScatter(N.filt, "S.Score", "nCount_RNA")
dev.off()

## CAN GET BETTER GENE LISTS, FOLLOW THIS TUTORIAL: https://hbctraining.github.io/scRNA-seq/lessons/05_SC_clustering_cells.html

## Running a PCA on cell cycle genes
all.genes <- rownames(N.filt)
N.filt <- ScaleData(N.filt, features = all.genes)
N.filt <- RunPCA(N.filt, features = c(s.genes, g2m.genes))
pdf("PCA_cellCycleGenes.pdf")
DimPlot(N.filt,reduction = "pca")
DimPlot(N.filt,reduction = "pca", group.by = 'trt')
dev.off()

table(paste0(N.filt$trt,N.filt$Phase))


N.filt=AddMetaData(N.filt, paste(N.filt$trt,N.filt$Phase, sep=";"),"group")

## Perform linear dimensional reduction
pdf("PCA_varGenes.pdf")
N.filt <- RunPCA(N.filt, features = VariableFeatures(object = N.filt), nfeatures.print = 5)
VizDimLoadings(N.filt, dims = 1:2, reduction = "pca")
DimPlot(N.filt, reduction = "pca", group.by = 'group', dims = c(1, 2),pt.size = 0.3)


DimPlot(N.filt, reduction = "pca", group.by = 'Phase')
DimHeatmap(N.filt, dims = 1:6, cells = 500, balanced = TRUE)
# Control that when running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase

## Determine the ‘dimensionality’ of the dataset - overcome technical noise
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
N.filt <- JackStraw(N.filt, num.replicate = 50, dims=30)
N.filt <- ScoreJackStraw(N.filt, dims = 1:30)
JackStrawPlot(N.filt, dims = 1:30)
## heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one
ElbowPlot(N.filt)
dev.off()

## Cluster the cells
N.filt <- FindNeighbors(N.filt, dims = 1:10)
N.filt <- FindClusters(N.filt, resolution = 0.5)  # resolution - reducing the number of clusters found in UMAP

## Run non-linear dimensional reduction (UMAP/tSNE)
N.filt <- RunUMAP(N.filt, reduction = "pca", dims = 1:6, n_neighbors = 10)

pdf("UMAPvarGenes.pdf")
DimPlot(N.filt, reduction = "umap", label=T)

DimPlot(N.filt, reduction = "umap", group.by = 'group', pt.size = 0.3)

## Highlight by Cell Cycle Phase
DimPlot(N.filt, reduction = "umap", group.by = 'Phase')

dev.off()

saveRDS(N.filt, file = "Nfilt.rds")
saveRDS(D, file = "N6dmso.rds")
saveRDS(W, file = "N6wee1i.rds")




