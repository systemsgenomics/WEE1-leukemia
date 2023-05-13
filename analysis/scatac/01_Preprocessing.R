######PREPROCESSING OF THE DATA#####

library(Signac); library(Seurat); library(GenomeInfoDb); library(EnsDb.Hsapiens.v75); library(ggplot2); library(hdf5r); library(data.table)

set.seed(1234)


##########T7#########
setwd('cellranger/T7/outs/')
## Seurat object setup

## Take only human peaks and name correctly for CreateGeneActivityMatrix function
# Use metadata to separate human and mouse cells, and make cell names unique
metadata <- read.csv(file = "singlecell.csv",
                     header = TRUE, row.names = 1)

id <- row.names(metadata)[as.logical(metadata$is_hg19_cell_barcode)&!(as.logical(metadata$is_mm10_cell_barcode))]  
hg <- metadata[rownames(metadata)%in%id,]

# added Wee1i to fragment file cell identifiers; same here
rownames(hg) <- paste0("Wee1i.", rownames(hg))

counts <- Read10X_h5(filename = ("filtered_peak_bc_matrix.h5"))
print(dim(counts))
counts <- counts[ ,colnames(counts)%in%id]
print(dim(counts))
colnames(counts) <- paste0("Wee1i.", colnames(counts))

## to avoid "-" prefix messing up coordinate extraction
## 
temp=rownames(counts)
rownames(counts)=sub("hg19_", "", rownames(counts))

mm10peaks=1:nrow(counts)%in%grep("mm10_", rownames(counts))
counts <- counts[!mm10peaks ,]
temp=temp[!mm10peaks]


Wee1i <- CreateSeuratObject(counts = counts, assay = "peaks", project = "Wee1i.atac", meta.data = hg)

rm(counts); rm(metadata)


## Sample T8
setwd('cellranger/T8/outs/')
metadata <- read.csv(file = "singlecell.csv",header = TRUE, row.names = 1)
id <- row.names(metadata)[as.logical(metadata$is_hg19_cell_barcode)&!(as.logical(metadata$is_mm10_cell_barcode))] 
hg <- metadata[rownames(metadata)%in%id,] 
# added Wee1i to fragment file cell identifiers; same here
rownames(hg) <- paste0("ctrl.", rownames(hg))

counts <- Read10X_h5(filename = ("filtered_peak_bc_matrix.h5"))
print(dim(counts))
counts <- counts[ ,colnames(counts)%in%id]
print(dim(counts))
colnames(counts) <- paste0("ctrl.", colnames(counts))

## hack to avoid "-" prefix messing up coordinate extraction
rownames(counts)=sub("hg19_", "", rownames(counts))


mm10peaks=1:nrow(counts)%in%grep("mm10_", rownames(counts))
counts <- counts[!mm10peaks ,]
#rownames(counts) <- gsub("mm10_","", rownames(counts)) #can be removed, but can also be run just in case

#rownames(counts) <- gsub("hg19_","", rownames(counts))

ctrl <- CreateSeuratObject(counts = counts, assay = "peaks", project = "ctrl.atac", meta.data = hg)
rm(counts); rm(metadata)
#head(colnames(Wee1i); head(rownames(Wee1i)); head(Wee1i@meta.data)[1:5];Wee1i

## ----------------- MERGE ----------------------------
## region-aware i.e. rename overlapping in ctrl using naming in Wee1
## https://satijalab.org/signac/reference/MergeWithRegions.html

# this will be our initial object but we will change it by adding feature quant from peaks
RS=MergeWithRegions(Wee1i, ctrl,sep.1 = c(":", "-"),sep.2 = c(":", "-"),distance=50)

write.table(rownames(Wee1i),file="analysis_T7_T8/temp_Wee1i_peaks.txt", sep="\t", row.names=F, col.names=F, quote=F)
write.table(rownames(ctrl),file="analysis_T7_T8/temp_ctrl_peaks.txt", sep="\t", row.names=F, col.names=F, quote=F)

## main script has the AWK + BEDTOOLS SCRIPT TO EDIT into BED and MERGE -> read in later in 01.2. preprocess to replace count matrix 

#head(colnames(Wee1i); head(rownames(Wee1i)); head(Wee1i@meta.data)[1:5];Wee1i

## leave out mouse cells from fragment file 
fragment.path <- 'analysis_T7_T8/T7.T8.fragments.tsv.gz'
fragment_file_filtered <- 'analysis_T7_T8/T7.T8.fragments_filt.tsv'
FilterFragments(fragment.path = fragment.path,cells = colnames(RS),output.path = fragment_file_filtered, assume.sorted=T )
# can check if worked with gunzip -c T7.T8.fragments_filt.tsv.bgz | cut -f 4 >tempids.txt

## set FRAGMENT INFO for Seurat object
fragment_file_filtered <- 'analysis_T7_T8/T7.T8.fragments_filt.tsv.bgz'
RS <- SetFragments(object = RS, file = fragment_file_filtered)


saveRDS(RS, file = 'analysis_T7_T8/T7_T8_prepro.rds')


