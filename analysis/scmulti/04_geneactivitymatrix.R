

#########CREATING A GENE ACTIVITY MATRIX#############

library(Signac); library(Seurat); library(GenomeInfoDb); library(EnsDb.Hsapiens.v86); library(data.table)



#To create a gene activity matrix, we extract gene coordinates for the human genome from EnsembleDB, and extend
#them to include the 2kb upstream region
RS <- readRDS('RS411_multiome_10h/RS_norm_dimred.rds')

gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding")

seqlevelsStyle(gene.coords) <- 'UCSC'

genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse', species="Homo sapiens")


genebodyandpromoter.coords <- Extend(x = genebody.coords, upstream = 2000, downstream = 0)

#head(genebodyandpromoter.coords)
class(genebodyandpromoter.coords)


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

RS[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities_RS)
DefaultAssay(RS) <- "ACTIVITY"
saveRDS(RS, file = 'RS411_multiome_10h/RS_genematrix.rds')



