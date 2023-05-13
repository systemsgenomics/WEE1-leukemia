
## conda activate R4_hdf5_seurat4
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(dittoSeq)

veh7="/MLL7vehicle.filt.rds"

#MLL7 treated
w1i7="/MLL7treated.filt.rds"

## load seurat objects
dd1=readRDS(veh7)
dd2=readRDS(w1i7)


mymeta=colnames(dd[[]])
myvars=c("Phase", "seurat_clusters","group")
myQC=c("nCount_RNA", "nFeature_RNA","percent.mt")

# PLOTS ----

minCol="slategrey"
maxCol="darkgoldenrod1"

## viz vars on umap
pdf("qcplots.pdf")
print(multi_dittoDimPlot(dd1, ncol=2,vars=myQC,  reduction.use="umap",min.col=minCol,max.col=maxCol,order="increasing",do.raster=T,raster.dpi=100,legend.size=3, legend.show=T,legend.title ="vehicle cells"))
print(multi_dittoDimPlot(dd2, ncol=2,vars=myQC,  reduction.use="umap",min.col=minCol,max.col=maxCol,order="increasing",do.raster=T,raster.dpi=100,legend.size=3, legend.show=T,legend.title ="treated cells"))
print(multi_dittoDimPlot(dd1, ncol=2,vars=myvars,  reduction.use="umap",min.col=minCol,max.col=maxCol,order="increasing",do.raster=T,raster.dpi=100,legend.size=3, legend.show=T,legend.title ="vehicle cells"))
print(multi_dittoDimPlot(dd2, ncol=2,vars=myvars,  reduction.use="umap",min.col=minCol,max.col=maxCol,order="increasing",do.raster=T,raster.dpi=100,legend.size=3, legend.show=T,legend.title ="treated cells"))
dev.off()


pdf("MLL7_ditto_dotplots2.pdf")
r.markers=c("CD34","MS4A1","CD19","CD22","CD72","CD81","PTEN","PPP3CC","PLCG2","RAC2","PIK3CD","POU2F2","ETS1","SPIB","MEF2D","EBF1","BCL6","NFKBID","SOCS1","SREBF1", "MYC")
dittoDotPlot(dd2, vars = r.markers, group.by = "seurat_clusters",scale = FALSE)
dev.off()
