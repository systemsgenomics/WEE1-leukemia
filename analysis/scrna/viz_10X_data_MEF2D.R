## --- viz features of interest from input data

## provide rds file as input
args <- commandArgs(TRUE)
rds=args[1]
name=sub(".rds", "", rds)
name=sub(".*\\/", "", name)

## conda activate R4_hdf5_seurat4
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(dittoSeq)

## ---
## function to perform dimensionality reduction and clustering
doDRclu=function(x,resol){
    #if(!isSCT) x <- ScaleData(x)
    #x <-RunPCA(x, ndims.print=1, nfeatures.print=1)
    ## k neighbors is here 20, UMAP done with k 15 and k 30
    x <- FindNeighbors(x, reduction = "pca", dims = 1:50,annoy.metric="cosine")
    x <- FindClusters(x, resolution = resol)
    ## version with preference to local structure preserving embedding
    x <- RunUMAP(x, reduction = "pca", dims = 1:50,n.neighbors=15, min.dist=0.3)
    x <- RunUMAP(x, reduction = "pca", dims = 1:50,n.neighbors=30, min.dist=0.5, reduction.name="uglobal")
    return(x)
}

## --

resDir="figures/"
if(!dir.exists(resDir)) dir.create(resDir)

## load seurat object
dd=readRDS(rds)
dd$scanpy_leidenr1_cluster=as.character(dd$scanpy_leidenr1_cluster)

s.genes <- c("MCM5","PCNA","TYMS", "FEN1","MCM2" ,"MCM4","RRM1","UNG", "GINS2","MCM6" ,"CDCA7" ,"DTL",
               "PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1", "GMNN" ,"WDR76" ,"SLBP","CCNE2",
               "UBR7","POLD3","MSH2", "ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM",
               "CASP8AP2","USP1","CLSPN","POLA1" , "CHAF1B","BRIP1","E2F8")
g2m.genes <- c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5",
                 "TPX2","TOP2A","NDC80", "CKS2", "NUF2","CKS1B","MKI67","TMPO","CENPF" ,"TACC3","FAM64A" ,"SMC4",
                 "CCNB2", "CKAP2L","CKAP2", "AURKB" ,"BUB1" ,"KIF11","ANP32E","TUBB4B" , "GTSE1" ,"KIF20B","HJURP","CDCA3",
                 "HN1" ,"CDC20","TTK" , "CDC25C","KIF2C","RANGAP1", "NCAPD2" ,"DLGAP5","CDCA2" , "CDCA8" ,"ECT2" ,"KIF23",
                 "HMMR", "AURKA","PSRC1","ANLN","LBR", "CKAP5" , "CENPE" , "CTCF" ,"NEK2" ,"G2E3" , "GAS2L3","CBX5","CENPA")

dd <- CellCycleScoring(object = dd, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)


## if this is asample w/o umap and cluster info then perform clustering and generate umap
if(!"umap"%in%getReductions(dd)) dd=doDRclu(dd,resol=0.5)


mymeta=colnames(dd[[]])
myvars=c("Phase", "seurat_clusters","scanpy_leidenr1_cluster")
myQC=c("nCount_RNA", "nFeature_RNA")


## define markers of interest
r.markers.w1=c("MEF2D","SREBF1","BCL6","ETS1", "SPIB","POU2F2","PTEN","PLCG2","PIK3CD")
r.markers.w1.2=c("IGLL1","VPREB1","CD79A","IGLL1", "VPREB3","CD79B","MEF2D","SREBF1","BCL6")
r.markers=list(r.markers.w1,r.markers.w1.2)
r.markers=r.markers[!lapply(r.markers,length)==0]



# PLOTS ----

minCol="slategrey"
maxCol="darkgoldenrod1"


## viz vars on umap for w1 paper
is_SCT=F
pdf("wee1_mef2d_dimplots.pdf"))
print(multi_dittoDimPlot(dd, ncol=2,vars=myvars, reduction.use="uglobal",legend.size=3, legend.show=T))
for(i in 1:length(r.markers)) print(multi_dittoDimPlot(dd, ncol=3,vars=r.markers[[i]], assay=ifelse(is_SCT,"SCT","RNA"), reduction.use="uglobal",min.col=minCol,max.col=maxCol,order="increasing",legend.size=3, legend.show=T))
dev.off()

