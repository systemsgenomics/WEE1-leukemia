##scATACseq samples - plotting from multiome atac

#module load r/3.6.3



library(Signac)
library(Seurat)

## motif scoring, AZD1775 treated cells that passed filtering

Wdm <- readRDS('Wee1icells.multi_072022_chromVAR_Homer.v2.rds')
DefaultAssay(Wdm) <- 'chromvar'


pdf("figures/10hATACumap_P53.Saos_motif_activity.pdf")
p1 <- FeaturePlot(Wdm,features = "p53.p53..Saos.p53.ChIP.Seq.GSE15780", reduction="umap",min.cutoff = 'q10',
max.cutoff = 'q90',pt.size = 0.5)
print(p1)
dev.off()

w1=readRDS("motifs_atac_lt_meta_RS411_wee1i.10h.rds")
metaw1=w1[[]]
lt1=metaw1$predicted.id
lt2=metaw1$predicted.idv2
lt=cbind(metaw1$barcode,lt1,lt2)
colnames(lt)=c("barcode", "RNAlt_predicted.id", "RNAlt_predicted.idv2")
lt=as.data.frame(lt,stringsAsFactors=F)

metaWdm=Wdm[[]]
tempLT=lt[match(metaWdm$barcode,lt$barcode),]
rownames(tempLT)=rownames(metaWdm)

Wdm=AddMetaData(Wdm,metadata=tempLT)
pdf("figures/10hATACumap_RNAlabelTransfer.pdf")
DimPlot(Wdm,group.by = 'RNAlt_predicted.id')
DimPlot(Wdm,group.by = 'RNAlt_predicted.idv2')
dev.off()

## -----


## motif scoring, DMSO treated cells
Ddm <- readRDS('DMSOcells.multi_072022_chromVAR_Homer.v2.rds')
## label transfer
ltD=read.table("barcode_to_predictedid_DMSO.txt", sep="\t", header=T, stringsAsFactors=F)


metaDdm=Ddm[[]]
tempLTD=ltD[match(metaDdm$barcode,ltD$barcode),]
rownames(tempLTD)=rownames(metaDdm)

Ddm=AddMetaData(Ddm,metadata=tempLTD)

pdf("figures/10hATACumap_RNAlabelTransfer_DMSO.pdf")
DimPlot(Ddm,group.by = 'RNAlt_predicted.id')
DimPlot(Ddm,group.by = 'RNAlt_predicted.idv2')
dev.off()

## -----


## fragment length plots

pdf('figures/fragLength_distr.pdf')

## all cells
f1=FragmentHistogram(
  RS,
  region = "chr1-1-2000000"
)
print(f1+ggtitle("all cells 10 h"))

## well matched PT11 cells
f2=FragmentHistogram(
  RS[,RS$isW==T&RS$RNAlt_predicted.idv2=="pt900"],
  region = "chr1-1-2000000"
)
print(f2+ggtitle("PT11 cells"))

## well matched PT9 cells
f3=FragmentHistogram(
  RS[,RS$isW==T&RS$RNAlt_predicted.idv2=="pt9"],
  region = "chr1-1-2000000"
)
print(f3+ggtitle("PT9 cells"))

## well matched PT5-7 cells
f4=FragmentHistogram(
  RS[,RS$isW==T&RS$RNAlt_predicted.idv2%in%c("pt5","pt6","pt7")],
  region = "chr1-1-2000000"
)
print(f4+ggtitle("PT5-7 cells"))



## DMSO cells
f5=FragmentHistogram(
  RS[,RS$isW==F],
  region = "chr1-1-2000000",
)
print(f5+ggtitle("DMSO cells"))

## other AZ clusters
f6=FragmentHistogram(
  RS[,RS$isW==T&RS$RNAlt_predicted.idv2!="pt900"&RS$RNAlt_predicted.idv2!="xNotReliable"],
  region = "chr1-1-2000000",
)
print(f6+ggtitle("AZ cells, not PT11"))

## DMSO cells
f7=FragmentHistogram(
  RS[,RS$isW==F&RS$RNAlt_predicted.idv2=="pt4"],
  region = "chr1-1-2000000",
)
print(f7+ggtitle("DMSO cells, PT4"))

f8=FragmentHistogram(
  RS[,RS$isW==F&RS$RNAlt_predicted.idv2%in%c("pt1","pt2","pt3")],
  region = "chr1-1-2000000",
)
print(f8+ggtitle("DMSO cells, PT1-3"))

f9=FragmentHistogram(
  RS[,RS$isW==F&RS$RNAlt_predicted.idv2=="pt900"],
  region = "chr1-1-2000000",
)
print(f9+ggtitle("DMSO cells, PT11"))


dev.off()



## atac fragments summaries 

## DMSO cells
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 #    85    6978   10258   11085   13254  257914

# PT5-7
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  106    7332   11027   12003   14964   64434

## PT9
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#  144    8918   12977   13931   17516   75335    1000

## PT11
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 163    6373    8295    9640   10548  106203    1000
