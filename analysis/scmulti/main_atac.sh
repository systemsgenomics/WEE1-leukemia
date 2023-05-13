## main script

## set input paths

PATH_W="M2/outs/"
PATH_D="M1/outs/"
WD="RS411_multiome_10h/"
## specify R version with Signac and Seurat installation; load it if running the R scripts
module load r/3.6.3
## load other needed tools, samtools contains bgzip and tabix needed for compressing/indexing fragment files
module load samtools/1.12
## bedtools
module load bedtools/2.27.1

date


## ---
## multiome ref genome has non-standard chromosomes: grep fragments corresponding to "proper" chromosomes
## needs to be done only once ; first check if completed
if [ ! -f $PATH_W"/chr_fragments.tsv" ]
then
cd $PATH_W
echo $PATH_W
gzip -cd atac_fragments.tsv.gz | grep chr | awk '{ if( $1 != "#" ){ print $0 } }' > chr_fragments.tsv
head chr_fragments.tsv

## add sample identifier into fragment data
awk 'BEGIN { OFS = "\t" } $14 = "Wee1i." $4' chr_fragments.tsv > temp.tsv
head temp.tsv

## reorganize columns and remove hg19_ prefix from chr names
awk -v OFS="\t" '{print $1, $2, $3, $6, $5}' temp.tsv  > chr_fragments.tsv
head chr_fragments.tsv
rm temp.tsv
fi

## -----------------------

## same for ctrl
if [ ! -f $PATH_D"/chr_fragments.tsv" ]
then
cd $PATH_D
echo $PATH_D
gzip -cd atac_fragments.tsv.gz | grep chr | awk '{ if( $1 != "#" ){ print $0 } }' > chr_fragments.tsv
head chr_fragments.tsv

awk 'BEGIN { OFS = "\t" } $14 = "ctrl." $4' chr_fragments.tsv > temp.tsv
##head temp.tsv

awk -v OFS="\t" '{print $1, $2, $3, $6, $5}' temp.tsv  > chr_fragments.tsv
head chr_fragments.tsv
rm temp.tsv
fi

if [ -f $PATH_D"/chr_fragments.tsv" ]
then
echo "fragment files with only fragments from standard chromosomes have been generated"
fi

## -----------------------
## merge fragments across treatments
## https://satijalab.org/signac/articles/faq.html?q=merge#how-do-i-merge-objects-with-signac
if [ ! -f $WD"/W.D.fragments.tsv.gz" ]
then
echo "merge fragments across treatments"



sort -k1,1V -k2,2n $PATH_W"/chr_fragments.tsv" $PATH_D"/chr_fragments.tsv" > $WD"/W.D.fragments.tsv"


# block gzip compress the merged file
# -@ 4 uses 4 threads
bgzip -@ 4 $WD"/W.D.fragments.tsv"

# index the bgzipped file
tabix --preset=bed $WD"/W.D.fragments.tsv.gz"
fi

if [ -f $WD"/W.D.fragments.tsv.gz" ]
then
echo "merged fragment file has been generated"
fi

## -----------------------

## process with Signac and Seurat

## preprocess samples (extract human cells, create Seurat objects)
R CMD BATCH 01_Preprocessing.R

## merge peaks from both samples so that can create a peak count matrix for the merged object
cat $PATH_W"/atac_peaks.bed" $PATH_D"/atac_peaks.bed" | sort -k1,1 -k2,2n > temp.sorted.bed
bedtools merge -i temp.sorted.bed > merged_peaks.bed
awk '{len = $3 - $2; if (len > 10) { print $_}}' merged_peaks.bed > merged_peaks_sizeFilt.bed

## region sizes before
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 192.0   848.0   895.0   869.5   918.0  2404.0
## after
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  192.0   848.0   895.0   869.5   918.0  2404.0

R CMD BATCH 01.2_Preprocessing.R

## calculate QC metrics and leave out poor quality cells
R CMD BATCH 02_qc-metrics.R

## normalize, cluster, find neighbors and plot the samples separately in lower dimension (UMAP)
R CMD BATCH 03_norm_dimred.R

## calculate gene activity matrix and save an object with this as default assay (ACTIVITY)
## the UMAP here is re-calculated using top genes; use this object to integrate with RNA-seq
R CMD BATCH 04_geneactivitymatrix.R


