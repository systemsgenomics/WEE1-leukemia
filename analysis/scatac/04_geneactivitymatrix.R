

#########CREATING A GENE ACTIVITY MATRIX#############

library(Signac); library(Seurat); library(GenomeInfoDb); library(EnsDb.Hsapiens.v75); library(data.table)



#To create a gene activity matrix, we extract gene coordinates for the human genome from EnsembleDB, and extend
#them to include the 2kb upstream region
RS <- readRDS('analysis_T7_T8/RS_norm_dimred.rds')

gene.coords <- genes(EnsDb.Hsapiens.v75, filter = ~ gene_biotype == "protein_coding")

seqlevelsStyle(gene.coords) <- 'UCSC'

genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse', species="Homo sapiens")


genebodyandpromoter.coords <- Extend(x = genebody.coords, upstream = 2000, downstream = 0)

head(genebodyandpromoter.coords)
class(genebodyandpromoter.coords)


#######T7#######
# create a gene by cell matrix


gene.activities_RS <- FeatureMatrix(
  
  fragments = GetFragments(object = RS, assay = 'peaks'),
  
  features = genebodyandpromoter.coords,
  
  cells = colnames(RS),
  
  chunk = 10
  
)


# convert rownames from chromsomal coordinates into gene names

gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)

rownames(gene.activities_RS) <- gene.key[rownames(gene.activities_RS)]

## filter genes without any signal
dim(gene.activities_RS)
## 20086
summary(rowSums(gene.activities_RS))
#gene.activities_wee1i = gene.activities_wee1i[rowSums(gene.activities_wee1i)>100,]
# this seems to drop very few genes so lets keep them for now

#dim(gene.activities_wee1i)
# add the gene activity matrix to the Seurat object as a new assay, and normalize it

RS[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities_RS)
DefaultAssay(RS) <- "ACTIVITY"
saveRDS(RS, file = 'analysis_T7_T8/RS_genematrix.rds')



RS <- RunTFIDF(RS, method=2)
RS <- FindTopFeatures(RS, min.cutoff = 'q0')
RS <- RunSVD(

object = RS,
n=150,

assay = 'ACTIVITY',

reduction.key = 'LSI.G_',

reduction.name = 'lsiG',
seed.use =1234

)

pdf("analysis_T7_T8/figures/genes_elbow_LSI.pdf")
ElbowPlot(RS, ndims = 150, reduction = "lsiG")
dev.off()

RS <- RunUMAP(object = RS, reduction = 'lsiG', dims = 1:50, seed.use = 1234)
RS@reductions$umap[[1:3, 1:2]]

RS <- FindNeighbors(object =RS, reduction = 'lsiG', dims = 1:50)
RS <- FindClusters(object = RS, verbose = FALSE, random.seed = 1234)


pdf('analysis_T7_T8/figures/clusters_geneactivity_RS.pdf')
DimPlot(object = RS, reduction = 'umap', label = TRUE) + NoLegend()
dev.off()



pdf('analysis_T7_T8/figures/ident_geneactivity.pdf')
DimPlot(RS, reduction = "umap", group.by = 'orig.ident')
dev.off()

pdf('analysis_T7_T8/figures/qc_geneactivity.pdf')
FeaturePlot(RS, features = c("nFeature_peaks", "TSS.enrichment", "blacklist_region_fragments","peak_region_fragments_mm10",
"enhancer_region_fragments", "promoter_region_fragments","GeneCountsPerCell"))
dev.off()



saveRDS(RS, file = 'analysis_T7_T8/RS_genematrix.rds')


