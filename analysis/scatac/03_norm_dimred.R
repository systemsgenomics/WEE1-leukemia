### object wiki https://github.com/satijalab/seurat/wiki/Seurat

#######NORMALIZATION AND DIMENSIONALITY REDUCTION#########

library(Seurat);library(Signac)

set.seed(1234)

RS <- readRDS('analysis_T7_T8/RS_qc.rds')

#######T7########

#Term frequency-inverse document frequency (TF-IDF) normalization of peaks by accessibility
RS <- RunTFIDF(RS, method=2)

## default method is 1
## 1: The LSI implementation used by Stuart & Butler et al. 2019 (https: //doi.org/10.1101/460147).
## 2: The standard LSI implementation used by Cusanovich & Hill et al. 2018 (https://doi.org/10.1016/j.cell.2018.06.052).
## 3: The log-TF method
## 4: The 10x Genomics method (no TF normalization)




#Feature selection, none applied since q0 keeps all
RS <- FindTopFeatures(RS, min.cutoff = 'q0')

#Dimensional reduction: Singular value decomposition (SVD) is run on the TD-IDF normalized matrix


RS <- RunSVD(

  object = RS,
  n=150,

  assay = 'peaks',

  reduction.key = 'LSI_',

  reduction.name = 'lsi',
  seed.use =1234

)



## examine result
#tt=RS@reductions$lsi
#tt[[1:3, 1:3]]
##LSI_1      LSI_2      LSI_3
##Wee1i.TTACTCAAGATCTAAG-1 2.042330 -0.5950781  0.3510433
##ctrl.AATGGCTTCGCTCGGA-1  4.630685 -3.9499145 -0.1060518
##ctrl.TACTGCCAGAGCAGCT-1  4.248084 -1.4000794  0.2221129


pdf("analysis_T7_T8/figures/elbow_LSI.pdf")
ElbowPlot(RS, ndims = 150, reduction = "lsi")
dev.off()


########## Non-linear dimension reduction and clustering ##############



RS <- RunUMAP(object = RS, reduction = 'lsi', dims = 1:50, seed.use = 1234)

##RS@reductions$umap[[1:3, 1:3]]
##UMAP_1    UMAP_2
##Wee1i.TTACTCAAGATCTAAG-1 -1.775942 -1.728919
##ctrl.AATGGCTTCGCTCGGA-1  -4.009875  1.019676
##ctrl.TACTGCCAGAGCAGCT-1   1.611668  2.912880


## k=20 by default
RS <- FindNeighbors(object =RS, reduction = 'lsi', dims = 1:50)


RS <- FindClusters(object = RS, verbose = FALSE, random.seed = 1234)
head(Idents(RS))
#Wee1i.TTACTCAAGATCTAAG-1  ctrl.AATGGCTTCGCTCGGA-1  ctrl.TACTGCCAGAGCAGCT-1
#7                        6                        4
#Wee1i.CAAGCTATCGTTACAG-1 Wee1i.GGTGCTGGTTCAGTTG-1 Wee1i.TAGCTTTGTAGAATAC-1
#5                        1                        1

pdf('analysis_T7_T8/figures/clusters_RS.260520.pdf')
DimPlot(object = RS, label = TRUE) + NoLegend()
## umap by default, then searches for tsne and then for pca
dev.off()


countsPerCell=colSums(as.matrix(GetAssayData(object = RS)))
print(summary(countsPerCell))
print(quantile(countsPerCell,0.9))
print(quantile(countsPerCell,0.1))
RS=AddMetaData(RS, metadata=countsPerCell, col.name="CountsPerCell")

pdf('analysis_T7_T8/figures/qcMap_RS.260520.pdf')
FeaturePlot(RS, features = c("nFeature_peaks", "TSS.enrichment", "blacklist_region_fragments","peak_region_fragments_mm10",
"enhancer_region_fragments", "promoter_region_fragments","CountsPerCell"))
dev.off()


## save original and subset
RS.old=RS

## compared to nuclei by CountsPerCell 0.1-0.9 quantile range
RS=subset(RS.old, subset = CountsPerCell > quantile(countsPerCell,0.1) & CountsPerCell < quantile(countsPerCell,0.9))
dim(RS)
## 168766   6600

RS.old=AddMetaData(RS.old, metadata=countsPerCell > quantile(countsPerCell,0.1) & countsPerCell < quantile(countsPerCell,0.9), col.name="subset")

## repeat SVD and downstream steps for subset
RS <- RunSVD(
object = RS,
assay = 'peaks',
reduction.key = 'LSI_',
reduction.name = 'lsi',
seed.use =1234
)

#tt=RS@reductions$lsi
#tt[[1:3, 1:3]]
##LSI_1      LSI_2       LSI_3
##ctrl.TACTGCCAGAGCAGCT-1  3.120576 -1.5802704 -0.04406963
##Wee1i.CAAGCTATCGTTACAG-1 2.393480  0.3435494  0.11559422
##Wee1i.GGTGCTGGTTCAGTTG-1 2.770689  5.1931683  0.01906964

########## Non-linear dimension reduction and clustering for the subset ##############



RS <- RunUMAP(object = RS, reduction = 'lsi', dims = 1:50, seed.use = 1234)
#RS@reductions$umap[[1:3, 1:2]]
##UMAP_1    UMAP_2
##ctrl.TACTGCCAGAGCAGCT-1   1.179552 -4.179198
##Wee1i.CAAGCTATCGTTACAG-1 -0.925856  2.563310
##Wee1i.GGTGCTGGTTCAGTTG-1 -4.338291  1.024973

RS <- FindNeighbors(object =RS, reduction = 'lsi', dims = 1:50)

RS <- FindClusters(object = RS, verbose = FALSE, random.seed = 1234)
head(Idents(RS))

pdf('analysis_T7_T8/figures/clusters_subset_RS.260520.pdf')

DimPlot(object = RS, reduction = 'umap', label = TRUE) + NoLegend()

dev.off()

pdf('analysis_T7_T8/figures/ident_subset_RS.260520.pdf')
DimPlot(RS, reduction = "umap", group.by = 'orig.ident')
dev.off()

pdf('analysis_T7_T8/figures/ident_RS.260520.pdf')
DimPlot(RS.old, reduction = "umap", group.by = 'orig.ident')
dev.off()

pdf('analysis_T7_T8/figures/subsetColored_RS.260520.pdf')
DimPlot(RS.old, reduction = "umap", group.by = 'subset')
dev.off()


pdf('analysis_T7_T8/figures/qc_subset_RS.260520.pdf')
FeaturePlot(RS, features = c("nFeature_peaks", "TSS.enrichment", "blacklist_region_fragments","peak_region_fragments_mm10",
"enhancer_region_fragments", "promoter_region_fragments","CountsPerCell"))
dev.off()


## --------------------------------
## save outputs

saveRDS(RS.old, file = 'analysis_T7_T8/RS_norm_dimred.rds')
saveRDS(RS, file = 'analysis_T7_T8/RS_norm_dimred_subset.rds')

RS2=RS
RS2 <- RunUMAP(object = RS2, reduction = 'lsi', dims = 2:50, seed.use = 1234)

pdf('analysis_T7_T8/figures/umap_woDim1.pdf')
DimPlot(RS2, reduction = "umap", group.by = 'orig.ident')
dev.off()

RS2 <- FindNeighbors(object =RS2, reduction = 'lsi', dims = 2:50)

RS2 <- FindClusters(object = RS2, verbose = FALSE, random.seed = 1234)

pdf('analysis_T7_T8/figures/umap_woDim1_clu.pdf')
DimPlot(RS2, reduction = "umap")
dev.off()
