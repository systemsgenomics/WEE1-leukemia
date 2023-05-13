#scRNAseq clusters comparison
#RS411 data cluster comparison DE analysis with

#https://www.rdocumentation.org/packages/Seurat/versions/3.1.1/topics/FindMarkers



library(Seurat)
library(dplyr)
library("ComplexHeatmap")

setwd("cluster_comp/")
RS.filt=readRDS("RSfilt.rds")


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


table(paste0(RS.filt$ptime,RS.filt$trt))

table(RS.filt$ptime)

pdf('ident_RS411_DMSO.wee1i_ptclusters.pdf')
DimPlot(RS.filt, reduction = "umap", group.by = 'ptime')
dev.off()

# cluster markers
markers <- FindAllMarkers(RS.filt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = 'RNA')
save(markers, file="markers.RData")
#load("markers.RData")

markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top1000 <- markers %>% group_by(cluster) %>% top_n(n = 1000, wt = avg_logFC)
class(top1000)
write.table(top1000, file="top1000clusterMarkers_RSfilt.txt", sep="\t", row.names=F, quote=F)

# find markers genes/statistics for scRNA RS411 wee1i treated cells
# comp.Find all markers for different wee1:n clusters. Compared to normal cell cycle G2M=PT4 (eg. clu3comp_clu5)

# find markers other way round cluster 1 specific 3
cluster1.markers3 <- FindMarkers(RS.filt, ident.1 = 1, ident.2 = 3, min.pct = 0.25)

# write cluster markers compared to cluster 1 to 3 out
write.table(cluster1.markers3, file="cluster1_comp_cluster3clusterMarkers.txt", sep="\t", row.names=T, quote=F)

# find markers other way round cluster 5 specific from 3
cluster5.markers3 <- FindMarkers(RS.filt, ident.1 = 5, ident.2 = 3, min.pct = 0.25)

# write cluster markers compared to cluster 5 to 3 out
write.table(cluster5.markers3, file="cluster5_comp_cluster3clusterMarkers.txt", sep="\t", row.names=T, quote=F)

# find markers other way round cluster 4 specific from 3
cluster4.markers3 <- FindMarkers(RS.filt, ident.1 = 4, ident.2 = 3, min.pct = 0.25)

# write cluster markers compared to cluster 4 to 3 out
write.table(cluster4.markers3, file="cluster4_comp_cluster3clusterMarkers.txt", sep="\t", row.names=T, quote=F)







#### RUN HEATMAP

#take mean expression by cluster pt time.
#first scale data

sdata=GetAssayData(RS.filt, slot="scale.data")
rownames(sdata)=rownames(RS.filt)

# enricher gene list from clu 3 (normal G2M) comp pt6 #
listbulk=unique(scan("cluster3_comp_cluster_all_wee1genes.txt",what="character"))
sellist=rownames(sdata)%in%listbulk

#take mean gene expression value inside that cluster. 
#take genes in correct order 
sdata=sdata[sellist,]
sdata=sdata[match(listbulk,rownames(sdata)),]
smeans=matrix(nrow=length(listbulk), ncol=11, data=0)
rownames(smeans)=listbulk
for(i in 1:nrow(smeans)){
  smeans[i,]=aggregate(sdata[i,]~RS.filt$ptime, FUN=mean)[,2]
}


#cluster rows using hiearchical clustering with 1-pearson correlation as distance
rowdendr=hclust(as.dist(1-cor(t(smeans),method="pearson")),method="average")

## do hiearchical clu yourself, pass to heatmap

## do kmeans yourself, pass to heatmap



set.seed(0)

pdf("clustercompheatmap.pdf", width = 5, height = 180) 
a=draw(Heatmap(smeans, km = 20, cluster_columns = F, column_names_gp = gpar(fontsize = 4),row_names_gp = gpar(fontsize = 4),column_title_gp = gpar(fontsize = 4)))
dev.off()
#print(a)



###
