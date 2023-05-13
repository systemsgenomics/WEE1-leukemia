###Integration of the scRNAseq data RS411 with Nalm-6 data ###

#https://satijalab.org/seurat/v3.0/atacseq_integration_vignette.html

library(Seurat); library(ggplot2)

setwd("/label_transf/")

## data sets

RS.filt <- readRDS("RSfilt.rds")
N.filt <- readRDS("Nfilt.rds")

RS.filt$tech <- "rna"
N.filt$tech <- "rna"

#pt-order cells

tempclu=as.vector(RS.filt[[]]$seurat_clusters)
tempclu[RS.filt[[]]$seurat_clusters==7]="pt1"
tempclu[RS.filt[[]]$seurat_clusters==2]="pt2"
tempclu[RS.filt[[]]$seurat_clusters==0]="pt3"
tempclu[RS.filt[[]]$seurat_clusters==3]="pt4"
tempclu[RS.filt[[]]$seurat_clusters==5]="pt5"
tempclu[RS.filt[[]]$seurat_clusters==4]="pt6"
tempclu[RS.filt[[]]$seurat_clusters==6]="pt7"
tempclu[RS.filt[[]]$seurat_clusters==10]="pt8"
tempclu[RS.filt[[]]$seurat_clusters==8]="pt9"
tempclu[RS.filt[[]]$seurat_clusters==9]="pt90"
tempclu[RS.filt[[]]$seurat_clusters==1]="pt900"

RS.filt[["ptime"]]=tempclu

colnames(RS.filt[[]])
colnames(N.filt[[]])



#Identifying the anchors between the scRNAseq RS411 ptime dataset and the scRNA-seq Nalm6 dataset and using these anchors to transfer the
#ptime labels from RS411 scRNA-seq data to the Nalm6 cells.

transfer.anchors <- FindTransferAnchors(reference = RS.filt, query = N.filt, features = VariableFeatures(object = RS.filt),
reference.assay = "RNA", query.assay = "RNA", reduction = "cca")

#To transfer the cluster ids, we provide a vector of previously annotated cell type labels for the RNA to the refdata parameter.
#The output will contain a matrix with predictions and confidence scores for each Nalm6 cells

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = RS.filt$ptime,weight.reduction = "cca")

N.filt <- AddMetaData(N.filt, metadata = celltype.predictions)

#We can then examine the distribution of prediction scores and optionally filter out those cells with low scores.

pdf('RNAintegration_prediction_scores_transfer_label_RS_PT_Nalm6.pdf')
hist(N.filt$prediction.score.max)
abline(v = 0.5, col = "red")
dev.off()

pdf('RNAintegration_predictedID_transfer_label_RS_PT_Nalm6.pdf')
DimPlot(N.filt, reduction = "umap", group.by = "predicted.id")
dev.off()

## create a separate annotation where poorly predicted cells are set as NotReliable
celltype.predictions2=celltype.predictions
colnames(celltype.predictions2)=paste(colnames(celltype.predictions2),2, sep="v")
celltype.predictions2$predicted.idv2[celltype.predictions2$prediction.score.maxv2 < 0.5]="xNotReliable"
N.filt <- AddMetaData(N.filt, metadata = celltype.predictions2)

pdf('RNAintegration_predictedID.onlyReliable_transfer_label_RS_PT_Nalm6.pdf')
DimPlot(N.filt, reduction = "umap", group.by = "predicted.idv2")
dev.off()


table(RS.filt$ptime)
table(N.filt$predicted.idv2)

write.table(N.filt[[]], file="meta.geneActivity.withPredictedRNAcluLabels_RS_PT_Nalm6_nocutoff.txt", sep="\t",quote=F)


## subset most reliable labels
table(N.filt$prediction.score.max > 0.2)
#USED 0.2 
N.filt.filtered <- subset(N.filt, subset = prediction.score.max > 0.2)


#####CO-EMBEDDING######

#this step is for visualization purposes only and is not a necessary part of the data transfer analysis.

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(RS.filt)
refdata <- GetAssayData(RS.filt, assay = "RNA", slot = "data")[genes.use, ]


# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the Nalm6 cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = "cca")

# this line adds the imputed data matrix to the Nalm6 data object
N.filt[["RNA"]] <- imputation


coembed <- merge(x = RS.filt, y = N.filt.filtered)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)


pdf('RNAmap_Nalm6RNAmap_integrated_RS_PT_Nalm6.pdf')
#CombinePlots(list(p1, p2))
DimPlot(coembed, reduction = "umap", group.by = 'tech')
DimPlot(coembed, reduction = "umap", group.by = 'ptime')
dev.off()

## can use this table to add labels to any other Nalm6 data object using match (not all cells are identical between objects)
write.table(N.filt[[]], file="meta.geneActivity.withPredictedRNAcluLabels_RS_PT_Nalm6.txt", sep="\t",quote=F)

saveRDS(coembed, file = "RNA_coembed_RS_PT_Nalm6.rds")

