#enhancers used for GROseq statistics (28k)
#module load r/3.6.3

library(chromVAR)

#########CREATING AN ENHANCER ACTIVITY MATRIX#############

library(Signac); library(Seurat); library(GenomeInfoDb); library(EnsDb.Hsapiens.v75); library(data.table)


## process coordinates for analysis
enhancers <- as.data.frame(read.table('enhancers_filt_centered.v3.bed6.bed',
header = F, sep = "\t", stringsAsFactors = F))
dim(enhancers)

colnames(enhancers)=c("chr", "start", "end", "seqname", "score", "strand")

## modify so that matches format of genebodyandpromoter.coords in script 04
## we need here the seqinfo object with chr and sizes
gene.coords <- genes(EnsDb.Hsapiens.v75, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse', species="Homo sapiens")

#chromosome name formatting
mychr=paste0("chr", 1:22)
mychr=c(mychr,"chrX", "chrY")
res1=keepSeqlevels(genebody.coords, mychr, pruning.mode = 'coarse')
res2=renameSeqlevels(res1, mychr)

seqinfo(test2)

enh=makeGRangesFromDataFrame(enhancers,
                         keep.extra.columns=FALSE,
                         ignore.strand=TRUE,
                         seqinfo=seqinfo(res2),
                         seqnames.field="chr",
                         start.field="start",
                         end.field=c("end"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)


# create an enhancers by cell matrix

RS <- readRDS('analysis_T7_T8/RS_qc.rds')


enhancers.activities_RS <- FeatureMatrix(

  fragments = GetFragments(object = RS, assay = 'peaks'),

  features = enh,

  cells = colnames(RS),

  chunk = 10

)


## check if enhancers without any signal
summary(rowSums(enhancers.activities_RS))
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##1.0   105.0   317.0   447.5   644.0  4920.0

RS[['EACTIVITY']] <- CreateAssayObject(counts = enhancers.activities_RS)
DefaultAssay(RS) <- "EACTIVITY"
saveRDS(RS, file = 'RS_enhancersmatrix.v2.rds')

countsPerCell=colSums(GetAssayData(object = RS, slot = "counts"))
print(summary(countsPerCell))
print(quantile(countsPerCell,0.9))
print(quantile(countsPerCell,0.1))

RS=AddMetaData(RS, metadata=countsPerCell, col.name="EnhancerCountsPerCell")
RS.old=RS
dim(RS)
## leave out cells that are having very low or very high counts
RS=subset(RS.old, subset = EnhancerCountsPerCell > quantile(countsPerCell,0.1) & EnhancerCountsPerCell < quantile(countsPerCell,0.9))
dim(RS)

RS <- NormalizeData(

  object = RS,

  assay = 'EACTIVITY',

  normalization.method = 'LogNormalize',

  scale.factor = median(colSums(GetAssayData(object = RS, slot = "counts")))

)



RS <- FindVariableFeatures(RS)
all.enh <- rownames(RS)
RS <- ScaleData(RS, features = all.enh)
RS <- RunPCA(RS,features = VariableFeatures(object = RS))
## Cluster the cells
RS <- FindNeighbors(RS, dims = 1:10)
RS <- FindClusters(RS, resolution = 0.5)  # resolution - reducing the number of clusters found in UMAP
## Run non-linear dimensional reduction (UMAP/tSNE)
RS <- RunUMAP(RS, reduction = "pca", dims = 1:6, n_neighbors = 10)


pdf('analysis_T7_T8/figures/clusters_enhancersactivity_RS.v2.pdf')
DimPlot(object = RS, reduction = 'umap', label = TRUE) + NoLegend()
dev.off()

pdf('analysis_T7_T8/figures/ident_enhanceractivity.v2.pdf')
DimPlot(RS, reduction = "umap", group.by = 'orig.ident')
dev.off()

pdf('analysis_T7_T8/figures/qc_enhanceractivity.v2.pdf')
FeaturePlot(RS, features = c("nFeature_peaks", "TSS.enrichment", "blacklist_region_fragments","peak_region_fragments_mm10",
"enhancer_region_fragments", "promoter_region_fragments","EnhancerCountsPerCell"))
dev.off()


saveRDS(RS, file = 'analysis_T7_T8/RS_enhancersmatrix.v2.rds')

## testing alternative dim red with LSI method
RS <- readRDS('analysis_T7_T8/RS_enhancersmatrix.v2.rds')

RS <- RunTFIDF(RS)

#Feature selection, none applied since q0 keeps all
RS <- FindTopFeatures(RS, min.cutoff = 'q0')

RS <- RunSVD(

object = RS,

assay = 'EACTIVITY',

reduction.key = 'LSI_',

reduction.name = 'lsi'

)

RS <- RunUMAP(object = RS, reduction = 'lsi', dims = 1:30)

RS <- FindNeighbors(object =RS, reduction = 'lsi', dims = 1:30)

RS <- FindClusters(object = RS, verbose = FALSE)



pdf('analysis_T7_T8/figures/enhancerclusters_LSI_RS.v2.pdf')

DimPlot(object = RS, label = TRUE) + NoLegend()
DimPlot(RS, reduction = "umap", group.by = 'orig.ident')
dev.off()
saveRDS(RS, file = 'RS_LSI_enhancersmatrix.v2.rds')
