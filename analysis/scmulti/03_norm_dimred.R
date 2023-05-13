### object wiki https://github.com/satijalab/seurat/wiki/Seurat

#######NORMALIZATION AND DIMENSIONALITY REDUCTION#########

library(Seurat);library(Signac)

set.seed(1234)

RS <- readRDS('RS411_multiome_10h/RS_qc.rds')

###############

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




pdf("RS411_multiome_10h/figures/elbow_LSI.pdf")
ElbowPlot(RS, ndims = 150, reduction = "lsi")
dev.off()


########## Non-linear dimension reduction and clustering ##############

## test leaving out component 1

RS <- RunUMAP(object = RS, reduction = 'lsi', dims = 1:50, seed.use = 1234)



## k=20 by default
RS <- FindNeighbors(object =RS, reduction = 'lsi', dims = 1:50)


saveRDS(RS, file = 'RS411_multiome_10h/RS_norm_dimred.rds')



pdf('RS411_multiome_10h/figures/clusters_RS.pdf')
DimPlot(object = RS, label = TRUE) + NoLegend()
## umap by default, then searches for tsne and then for pca
dev.off()



pdf('RS411_multiome_10h/qcMap_RS.pdf')
FeaturePlot(RS, features = c("atac_peak_region_fragments", "TSS.enrichment", "atac_TSS_fragments","nucleosome_signal","pct_reads_in_peaks"))
dev.off()

