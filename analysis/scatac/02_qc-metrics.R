
library(Signac); library(Seurat); library(GenomeInfoDb); library(EnsDb.Hsapiens.v75); library(ggplot2); library(data.table)


RS <- readRDS('analysis_T7_T8/RS_prepro.rds')
dim(RS)


######COMPUTING QC-METRICS########


# create granges object with TSS positions

gene.ranges <- genes(EnsDb.Hsapiens.v75)
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]

tss.ranges <- GRanges(
seqnames = seqnames(gene.ranges),
ranges = IRanges(start = start(gene.ranges), width = 2),
strand = strand(gene.ranges)
)

seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')
## rename chr with hg19_ prefix

#### QC #####


RS <- NucleosomeSignal(object = RS, region='chr1-1-249250621')



RS$pct_reads_in_peaks <- RS$peak_region_fragments /RS$total * 100

RS$blacklist_ratio <- RS$blacklist_region_fragments / RS$peak_region_fragments



plot1 <- VlnPlot(
  
  object = RS,
  
  features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal'),
  
  pt.size = 0.1) + NoLegend()



plot2_a <- VlnPlot(
  
  object = RS,
  
  features = 'peak_region_fragments',
  
  pt.size = 0.1, log = TRUE) + NoLegend()



plot2_b <- FeatureScatter(RS,"peak_region_fragments",'nucleosome_signal', pt.size = 0.1) + NoLegend()

plot2_c <- FeatureScatter(RS,"peak_region_fragments",'blacklist_ratio', pt.size = 0.1) + NoLegend()

plot2 <- CombinePlots(plots = list(plot2_a,plot2_b,plot2_c), ncol = 3)



pdf('analysis_T7_T8/figures/qc_RS.pdf')

CombinePlots(list(plot1,plot2),ncol = 1)

dev.off()


##group by cells with high or low nucleosomal signal strength.


RS$nucleosome_group <- ifelse(RS$nucleosome_signal > 10, 'NS > 10', 'NS < 10')

#PeriodPlot(object = Wee1i, group.by = 'nucleosome_group', region = "hg19_chr1-1-2000000")

pdf('analysis_T7_T8/figures/ns10_RS.pdf')

PeriodPlot(object =RS, group.by = 'nucleosome_group', region = "chr1-1-2000000")

dev.off()



# to save time use the first 2000 TSSs
RS <- TSSEnrichment(object = RS, tss.positions = tss.ranges[1:2000])

RS$high.tss <- ifelse(RS$TSS.enrichment > 2, 'High', 'Low')

pdf('figures/high.tss_RS.pdf')

TSSPlot(RS, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()

dev.off()

isW=1:nrow(RS[[]])%in%grep("Wee1i",rownames(RS[[]]))
#FALSE  TRUE
# 8053 10206
RS=AddMetaData(RS,isW, col.name="isW")


write.table(RS[[]], file="RS_preqc_meta.txt", sep="\t", row.names=F, col.names=T, quote=F)


## total cells 4446
#filter out cells that are outliers for the QC metrics.
#length(which(Wee1i$peak_region_fragments>1000)) #4177
#length(which(Wee1i$peak_region_fragments< 20000)) #2184
## CHANGED TO length(which(Wee1i$peak_region_fragments< 100000)) #4440
#length(which(Wee1i$pct_reads_in_peaks>15)) #3930
#length(which(Wee1i$blacklist_ratio<0.05)) #4446
#length(which(Wee1i$nucleosome_signal<10)) #4093
#length(which(Wee1i$TSS.enrichment>2)) #3588
dim(RS)

## subset to get rid of most likely artefacts but w/o using metrics that differ in mitotic cell subset (PT11 matched cell features based on 10h multiome result)
RS2 <- subset(RS, subset = peak_region_fragments > 1000 & peak_region_fragments < 100000 & blacklist_ratio < 0.05 &  TSS.enrichment > 2& nucleosome_signal < 50)

pdf("figures/split_by_sample_vio.pdf")
VlnPlot( object = RS2, features = c('pct_reads_in_peaks',  'nucleosome_signal'),
  pt.size = 0.1, split.by="isW") + NoLegend()
VlnPlot( object = RS2[,RS2$isW==T], features = c('nucleosome_signal'),pt.size = 0)
VlnPlot( object = RS2[,RS2$isW==F], features = c('nucleosome_signal'),pt.size = 0)
dev.off()


pdf("figures/split_by_sample_vio.v2.pdf")
VlnPlot( object = RS2, features = c('pct_reads_in_peaks',  'nucleosome_signal'),
  pt.size = 0, split.by="isW") + NoLegend()
VlnPlot( object = RS2[,RS2$isW==T], features = c('pct_reads_in_peaks',  'nucleosome_signal'),
  pt.size = 0) + NoLegend()
VlnPlot( object = RS2[,RS2$isW==F], features = c('pct_reads_in_peaks',  'nucleosome_signal'),
  pt.size = 0) + NoLegend()
dev.off()


RS <- subset(RS, subset = peak_region_fragments > 1000 & peak_region_fragments < 100000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
dim(RS)
#3260 cells remain





## save objects
saveRDS(RS, file = 'analysis_T7_T8/RS_qc.rds')



