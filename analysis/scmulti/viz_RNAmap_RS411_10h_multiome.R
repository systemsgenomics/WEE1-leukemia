##scRNAseq samples processed separately - plotting features from multiome atac

#conda activate R4_hdf5_seurat4
#
library(Seurat)
library(dplyr)


## labe transfer output object, contains the 10 h sample and coembedded (imputed) values for 24 h reference
w1meta=read.table("meta.RS41124h.geneActivity.withPredictedRNAcluLabels_RS411_10h.wee1i.txt", sep="\t", header=T, stringsAsFactors=F)
dmeta=read.table("meta.RS41124h.geneActivity.withPredictedRNAcluLabels_RS411_10h.DMSO.txt", sep="\t", header=T, stringsAsFactors=F)

setwd("scMultiome/")
w1=readRDS("RS411_wee1i.10h.rds")
d=readRDS("RS411_DMSO.10h.rds")

## wo cutoffs
w1 <- AddMetaData(w1, metadata = w1meta$predicted.id,col.name="predicted.id")
d <- AddMetaData(d, metadata = dmeta$predicted.id,col.name="predicted.id")

w1 <- AddMetaData(w1, metadata = w1meta$predicted.idv2,col.name="predicted.idv2")
d <- AddMetaData(d, metadata = dmeta$predicted.idv2,col.name="predicted.idv2")


setwd("RS411_multiome_10h/")

## read in umap coords generated with scanpy (velocyto run)
uw=read.table("M2_umap.csv", sep=",", header=T, stringsAsFactors=F)
ud=read.table("M1_umap.csv", sep=",", header=T, stringsAsFactors=F)

all(sub("_.*","",uw[,1])==sub("-.*","",colnames(w1)))
all(sub("_.*","",ud[,1])==sub("-.*","",colnames(d)))

rownames(uw)=sub("_.*","-1_2",uw[,1])
rownames(ud)=sub("_.*","-1_1",ud[,1])

uw=uw[,-1]
ud=ud[,-1]

colnames(ud)=c("UMAP_1", "UMAP_2")
colnames(uw)=c("UMAP_1", "UMAP_2")

# store this as a custom dimensional reduction
w1[["scanpyumap"]] <-  CreateDimReducObject(
  embeddings = as.matrix(uw),
  loadings = matrix(data=1,nrow=dim(w1)[1],ncol=2),
  stdev = rep(1,dim(w1)[1]),
  key = "scanpyumap_",
  assay = "RNA"
)

d[["scanpyumap"]] <-  CreateDimReducObject(
  embeddings = as.matrix(ud),
  loadings = matrix(data=1,nrow=dim(d)[1],ncol=2),
  stdev = rep(1,dim(d)[1]),
  key = "scanpyumap_",
  assay = "RNA"
)

pdf("figures/RNAmaps_lt24h.pdf")
DimPlot(w1, reduction = "umap", group.by = 'predicted.id')
DimPlot(d, reduction = "umap", group.by = 'predicted.id')
DimPlot(w1, reduction = "umap", group.by = 'predicted.idv2')
DimPlot(d, reduction = "umap", group.by = 'predicted.idv2')
dev.off()

pdf("figures/RNAmaps_lt24h_scanpyumap.pdf")
DimPlot(w1, reduction = "scanpyumap", group.by = 'predicted.id')
DimPlot(d, reduction = "scanpyumap", group.by = 'predicted.id')
DimPlot(w1, reduction = "scanpyumap", group.by = 'predicted.idv2')
DimPlot(d, reduction = "scanpyumap", group.by = 'predicted.idv2')
dev.off()


pdf("figures/RNAmaps_cellcycle_scanpyumap.pdf")
DimPlot(w1, reduction = "scanpyumap", group.by = 'Phase')
dev.off()

pdf("figures/RNAmaps_Sscore_scanpyumap.pdf")
FeaturePlot(w1, reduction = "scanpyumap", features = 'S.Score')
dev.off()

pdf("figures/RNAmaps_G2Mscore_scanpyumap.pdf")
FeaturePlot(w1, reduction = "scanpyumap", features = 'G2M.Score')
dev.off()



w1 <- AddMetaData(w1, metadata = sub("_2", "", colnames(w1)),col.name="barcode")
d <- AddMetaData(d, metadata = sub("_1", "", colnames(d)),col.name="barcode")

## multiome cellranger outputs
metaW=read.table("M2/outs/per_barcode_metrics.csv", sep=",", header=T, stringsAsFactors=F)

metaD=read.table("M1/outs/per_barcode_metrics.csv", sep=",", header=T, stringsAsFactors=F)

metaW=metaW[metaW$barcode%in%w1$barcode,]
metaD=metaD[metaD$barcode%in%d$barcode,]

metaW=metaW[match(w1$barcode,metaW$barcode),]
metaD=metaD[match(d$barcode,metaD$barcode),]

## additional metadata from RS411_multiome_10h/RS_preqc.rds
qcA=read.table("RS_preqc_meta.txt", sep="\t",header=T, stringsAsFactors=F)
qcW=qcA[qcA$isW==T,]
qcD=qcA[qcA$isW==F,]

tempW=qcW[match(metaW$barcode,qcW$barcode),]
tempD=qcD[match(metaD$barcode,qcD$barcode),]

metaW$pct_reads_in_peaks=tempW$pct_reads_in_peaks
metaD$pct_reads_in_peaks=tempD$pct_reads_in_peaks

metaW$nucleosome_signal=tempW$nucleosome_signal
metaD$nucleosome_signal=tempD$nucleosome_signal

metaW$TSS.enrichment=tempW$TSS.enrichment
metaD$TSS.enrichment=tempD$TSS.enrichment

rownames(metaW)=rownames(w1[[]])
rownames(metaD)=rownames(d[[]])

w1 <- AddMetaData(w1, metadata = metaW)
d <- AddMetaData(d, metadata = metaD)

pdf("figures/RNAmaps_atacQC.pdf")
FeaturePlot(w1, reduction = "umap", features = c('pct_reads_in_peaks','nucleosome_signal','TSS.enrichment'))
FeaturePlot(d, reduction = "umap", features = c('pct_reads_in_peaks','nucleosome_signal','TSS.enrichment'))
dev.off()
