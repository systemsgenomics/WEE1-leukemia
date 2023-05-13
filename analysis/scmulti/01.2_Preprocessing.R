library(Signac); library(Seurat); library(GenomeInfoDb); library(EnsDb.Hsapiens.v75); library(ggplot2); library(data.table); library(Matrix)

#' BED to GRanges
#'
#' This function loads a BED-like file and stores it as a GRanges object.
#' The tab-delimited file must be ordered as 'chr', 'start', 'end', 'id', 'score', 'strand'.
#' The minimal BED file must have the 'chr', 'start', 'end' columns.
#' Any columns after the strand column are ignored.
#'
#' @param file Location of your file
#' @keywords BED GRanges
#' @export
#' @examples
#' bed_to_granges('my_bed_file.bed')

bed_to_granges <- function(file){
    df <- read.table(file,
    header=F,
    stringsAsFactors=F, sep="\t")
    
    if(length(df) > 6){
        df <- df[,-c(7:length(df))]
    }
    
    if(length(df)<3){
        stop("File has less than 3 columns")
    }
    
    header <- c('chr','start','end','id','score','strand')
    names(df) <- header[1:length(names(df))]
    
    if('strand' %in% colnames(df)){
        df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
    }
    
    
    
    if(length(df)==3){
        gr <- with(df, GRanges(chr, IRanges(start, end)))
    } else if (length(df)==4){
        gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
    } else if (length(df)==5){
        gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
    } else if (length(df)==6){
        gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
    }
    return(gr)
}


## get merged peaks and convert to Granges
mpeaks=bed_to_granges('RS411_multiome_10h/merged_peaks_sizeFilt.bed')

## object with cells to keep
#RS <- readRDS('RS411_multiome_10h/W_D_prepro.rds')

### re-calculate count matrix for merged object

peaks_merged <- FeatureMatrix(
fragments = 'RS411_multiome_10h/W.D.fragments.tsv.gz',
features = mpeaks,
chunk = 10
)
#cells = colnames(RS)

print(ncol(peaks_merged))
head(colnames(peaks_merged))

save(peaks_merged, file="temp.RData")

## filtered file contains more cells, subset to those in metadata marked as cells

metaD=read.table("M1/outs/per_barcode_metrics.csv", sep=",", header=T, stringsAsFactors=F)
metaD=metaD[metaD$is_cell==1,]

head(metaD$atac_barcode)
rownames(metaD)=paste("ctrl.",metaD$barcode, sep="")
head(rownames(metaD))


metaW=read.table("M2/outs/per_barcode_metrics.csv", sep=",", header=T, stringsAsFactors=F)
metaW=metaW[metaW$is_cell==1,]
rownames(metaW)=paste("Wee1i.",metaW$barcode, sep="")

meta=rbind(metaD, metaW)

rs_peaks=peaks_merged[,which(colnames(peaks_merged)%in%rownames(meta))]

print(ncol(rs_peaks))
head(rownames(rs_peaks))

meta2=meta[match(colnames(rs_peaks),rownames(meta)),]
all(colnames(rs_peaks)==rownames(meta2))


num_cells_ncounted =tabulate(rs_peaks@i + 1)

## peak needs to have count in at least 1 % of cells
th=ncol(rs_peaks)*0.01
dim(rs_peaks)
rs_peaks=rs_peaks[which(num_cells_ncounted>th),]
## from 210 k peaks to 95 k peaks
dim(rs_peaks)

RS2 <- CreateSeuratObject(counts = rs_peaks, assay = "peaks", project = "RS.atac", meta.data=meta2)
RS2 <- SetFragments(object = RS2, file = 'RS411_multiome_10h/W.D.fragments.tsv.gz')


saveRDS(RS2, file = 'RS411_multiome_10h/RS_prepro.rds')


