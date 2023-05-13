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
mpeaks=bed_to_granges('analysis_T7_T8/merged_peaks_sizeFilt.bed')

## object with cells to keep
RS <- readRDS('analysis_T7_T8/T7_T8_prepro.rds')

### re-calculate count matrix for merged object

peaks_merged <- FeatureMatrix(
fragments = 'analysis_T7_T8/T7.T8.fragments_filt.tsv.bgz',
features = mpeaks,
chunk = 10
)
#cells = colnames(RS)

print(ncol(peaks_merged))
## just in case the filtered file contains more cells, subset again
rs_peaks=peaks_merged[,which(colnames(peaks_merged)%in%colnames(RS))]
print(ncol(rs_peaks))
head(rownames(rs_peaks))

binary_mat = as.matrix((rs_peaks > 0) + 0)
binary_mat = Matrix(binary_mat, sparse = TRUE)
binary_mat[1:3,1:3]
num_cells_ncounted = rowSums(binary_mat)

## can save here if relevant for other analyses
#saveRDS(binary_mat[which(num_cells_ncounted>ncol(binary_mat)*0.01),],file='analysis_T7_T8/binary_mat.rds')

## peak needs to have count in at least 1 % of cells
rs_peaks=rs_peaks[which(num_cells_ncounted>ncol(binary_mat)*0.01),]

RS2 <- CreateSeuratObject(counts = rs_peaks, assay = "peaks", project = "RS.atac", meta.data = RS[[]])
RS2 <- SetFragments(object = RS2, file = 'analysis_T7_T8/T7.T8.fragments_filt.tsv.bgz')


saveRDS(RS2, file = 'analysis_T7_T8/RS_prepro.rds')

