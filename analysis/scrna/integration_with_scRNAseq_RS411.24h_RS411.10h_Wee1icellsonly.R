###Integration of the scRNAseq data with multiome RS411 10h wee1i cells object###


library(Seurat); library(ggplot2)

setwd("label.transfer")

## data sets

#24h
RS.filt <- readRDS("RSfilt.rds")

#RS411 10h multiome wee1i. cells only GEX (first run ) called here N.filt
N.filt <- readRDS("RS411_wee1i.10h.rds")

RS.filt$tech <- "rna"
N.filt$tech <- "rna"

# pt-order cells

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




#Identifying the anchors between the scRNAseq RS411 ptime dataset and the scRNA-seq RS 10 h dataset and using these anchors to transfer the
#ptime labels from RS scRNA-seq data to the Nalm6 cells.

transfer.anchors <- FindTransferAnchors(reference = RS.filt, query = N.filt, features = VariableFeatures(object = RS.filt),
reference.assay = "RNA", query.assay = "RNA", reduction = "cca")





#To transfer the cluster ids, we provide a vector of previously annotated cell type labels for the RNA to the refdata parameter.
#The output will contain a matrix with predictions and confidence scores for each Nalm6 cells

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = RS.filt$ptime,weight.reduction = "cca")


N.filt <- AddMetaData(N.filt, metadata = celltype.predictions)

#We can then examine the distribution of prediction scores and optionally filter out those cells with low scores.

pdf('RNAintegration_prediction_scores_transfer_label.wee1i.pdf')
hist(N.filt$prediction.score.max)
abline(v = 0.5, col = "red")
dev.off()

pdf('RNAintegration_predictedID_transfer_label.wee1i.pdf')
DimPlot(N.filt, reduction = "umap", group.by = "predicted.id")
dev.off()

## create a separate annotation where poorly predicted cells are set as NotReliable
celltype.predictions2=celltype.predictions
colnames(celltype.predictions2)=paste(colnames(celltype.predictions2),2, sep="v")
celltype.predictions2$predicted.idv2[celltype.predictions2$prediction.score.maxv2 < 0.5]="xNotReliable"
N.filt <- AddMetaData(N.filt, metadata = celltype.predictions2)

pdf('RNAintegration_predictedID.onlyReliable_transfer_label.wee1i.pdf')
DimPlot(N.filt, reduction = "umap", group.by = "predicted.idv2")
dev.off()

table(celltype.predictions2$predicted.idv2[celltype.predictions2$prediction.score.maxv2 < 0.5])


table(N.filt$predicted.idv2)



table(RS.filt$ptime)


table(N.filt$trt)


table(N.filt$prediction.score.max > 0.2)

#with > 0.2 keep all the cells. label transfer went well

N.filt.filtered <- subset(N.filt, subset = prediction.score.max > 0.2)


#####CO-EMBEDDING######

#this step is for visualization purposes only and is not a necessary part of the data transfer analysis.

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(RS.filt)
refdata <- GetAssayData(RS.filt, assay = "RNA", slot = "data")[genes.use, ]


# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the 10h data cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = "cca")

# this line adds the imputed data matrix to the 10h data object
N.filt[["RNA"]] <- imputation

#0.2 cutoff used for RS411 10h

coembed <- merge(x = RS.filt, y = N.filt.filtered)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

###what celltype is replaced with here? --> seurat_clusters?
#coembed$ptime <- ifelse(!is.na(coembed$ptime), coembed$ptime)



pdf('RNAmap_PTlabels_integrated.wee1i.pdf')
#CombinePlots(list(p1, p2))
DimPlot(coembed, reduction = "umap", group.by = 'tech')
DimPlot(coembed, reduction = "umap", group.by = 'ptime')
dev.off()

## can use this table to add labels to any other RS411 10h data  object using match (not all cells are identical between objects)
write.table(N.filt[[]], file="meta.RS41124h.geneActivity.withPredictedRNAcluLabels_RS411_10h.wee1i.txt", sep="\t",quote=F)

saveRDS(coembed, file = "RNA_coembed.wee1i.rds")


####
