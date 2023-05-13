
library(Signac); library(Seurat); library(GenomeInfoDb); library(EnsDb.Hsapiens.v86); library(ggplot2); library(data.table)


RS <- readRDS('RS411_multiome_10h/RS_prepro.rds')
dim(RS)


######COMPUTING QC-METRICS########


# create granges object with TSS positions

gene.ranges <- genes(EnsDb.Hsapiens.v86)
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

#Calculate the strength of the nucleosome signal per cell. Computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) to fragments < 147 bp (nucleosome-free)
RS <- NucleosomeSignal(object = RS, region='chr1-1-249250621')



RS$pct_reads_in_peaks <- RS$atac_peak_region_fragments /RS$atac_fragments * 100

#RS$blacklist_ratio <- RS$blacklist_region_fragments / RS$peak_region_fragments

#"atac_raw_reads"
#[23] "atac_unmapped_reads"               "atac_lowmapq"
#[25] "atac_dup_reads"                    "atac_chimeric_reads"
#[27] "atac_mitochondrial_reads"          "atac_fragments"
#[29] "atac_TSS_fragments"                "atac_peak_region_fragments"
#[31] "atac_peak_region_cutsites"


plot1 <- VlnPlot(
  object = RS,
  features = c('pct_reads_in_peaks', 'atac_dup_reads', 'nucleosome_signal'),
  pt.size = 0.1) + NoLegend()



plot2_a <- VlnPlot(
  object = RS,
  features = 'atac_peak_region_fragments',
  pt.size = 0.1, log = TRUE) + NoLegend()



plot2_b <- FeatureScatter(RS,"atac_peak_region_fragments",'nucleosome_signal', pt.size = 0.1) + NoLegend()
#plot2_c <- FeatureScatter(RS,"atac_peak_region_fragments",'blacklist_ratio', pt.size = 0.1) + NoLegend()
plot2 <- CombinePlots(plots = list(plot2_a,plot2_b), ncol = 2)



pdf('RS411_multiome_10h/figures/qc_RS.pdf')
CombinePlots(list(plot1,plot2),ncol = 1)
dev.off()


##group by cells with high or low nucleosomal signal strength.
RS$nucleosome_group <- ifelse(RS$nucleosome_signal > 10, 'NS > 10', 'NS < 10')

pdf('RS411_multiome_10h/figures/ns10_RS.pdf')
PeriodPlot(object =RS, group.by = 'nucleosome_group', region = "chr1-1-2000000")
dev.off()



# to save time use the first 2000 TSSs
RS <- TSSEnrichment(object = RS, tss.positions = tss.ranges[1:2000])
RS$high.tss <- ifelse(RS$TSS.enrichment > 2, 'High', 'Low')

pdf('RS411_multiome_10h/figures/high.tss_RS.pdf')
TSSPlot(RS, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.off()

RS$high.tss5 <- ifelse(RS$TSS.enrichment > 5, 'High', 'Low')

pdf('RS411_multiome_10h/figures/high.tss5_RS.pdf')
TSSPlot(RS, group.by = 'high.tss5') + ggtitle("TSS enrichment score") + NoLegend()
dev.off()



dim(RS)
## not using & blacklist_ratio < 0.05, changed pct reads in peaks min 50
## subset by QC metrics
RS <- subset(RS, subset = atac_peak_region_fragments > 1000 & atac_peak_region_fragments < 100000 & pct_reads_in_peaks > 50 & nucleosome_signal < 10 )
dim(RS)
# 24297


## save objects
saveRDS(RS, file = 'RS411_multiome_10h/RS_qc.rds')

write.table(RS[[]], file="RS411_multiome_10h/qcFilt_atac_metadata.txt",sep="\t", row.names=F, col.names=F, quote=F)


## COMPARE HERE TO gexp data passing filters:
#GEX=readRDS("scMultiome/RSfilt.rds")
#gID=sub("_1","",rownames(GEX[[]]))
#gID=sub("_2","",rownames(GEX[[]]))
#length(which(gID%in%RS$barcode))
##[1] 11403
