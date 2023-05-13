######PREPROCESSING OF THE DATA#####
#module load r/3.6.3
library(Signac); library(Seurat); library(GenomeInfoDb); library(EnsDb.Hsapiens.v75); library(ggplot2); library(hdf5r); library(data.table)

set.seed(1234)


##########10 h WEE1i RS4;11 #########
setwd("M2/outs/")
## Seurat object setup
counts <- Read10X_h5(filename = ("filtered_peak_bc_matrix.h5"))
print(dim(counts))
colnames(counts) <- paste0("Wee1i.", colnames(counts))

Wee1i <- CreateSeuratObject(counts = counts, assay = "peaks", project = "Wee1i.atac")
rm(counts)


## DMSO RS4;11 10 h
setwd('M1/outs/')
counts <- Read10X_h5(filename = ("filtered_peak_bc_matrix.h5"))
print(dim(counts))
colnames(counts) <- paste0("ctrl.", colnames(counts))

ctrl <- CreateSeuratObject(counts = counts, assay = "peaks", project = "ctrl.atac")
rm(counts)

## ----------------- MERGE ----------------------------
## region-aware i.e. rename overlapping in ctrl using naming in Wee1
## https://satijalab.org/signac/reference/MergeWithRegions.html

# this will be our initial object but we will change it by adding feature quant from peaks
RS=MergeWithRegions(Wee1i, ctrl,sep.1 = c(":", "-"),sep.2 = c(":", "-"),distance=50)

write.table(rownames(Wee1i),file="temp_Wee1i_peaks.txt", sep="\t", row.names=F, col.names=F, quote=F)
write.table(rownames(ctrl),file="temp_ctrl_peaks.txt", sep="\t", row.names=F, col.names=F, quote=F)

## main script has the AWK + BEDTOOLS SCRIPT TO EDIT into BED and MERGE -> read in later in 01.2. preprocess to replace count matrix 


fragment.path <- 'RS411_multiome_10h/W.D.fragments.tsv.gz'
fragment_file_filtered <- 'W.D.fragments_filt.tsv'
FilterFragments(fragment.path = fragment.path,cells = colnames(RS),output.path = fragment_file_filtered, assume.sorted=T )

## set FRAGMENT INFO for Seurat object
fragment_file_filtered <- 'RS411_multiome_10h/W.D.fragments_filt.tsv.bgz'
RS <- SetFragments(object = RS, file = fragment_file_filtered)


saveRDS(RS, file = 'RS411_multiome_10h/W_D_prepro.rds')

