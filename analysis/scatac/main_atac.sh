## main script

## set input paths

PATH_T7="cellranger/T7/outs/"
PATH_T8="cellranger/T8/outs/"
WD="analysis_T7_T8/"
## specify R version with Signac and Seurat installation; load it if running the R scripts
#module load r/3.6.0
## load other needed tools, samtools contains bgzip and tabix needed for compressing/indexing fragment files
#module load samtools/1.9
## bedtools
#module load bedtools2/2.27.1

date

## mouse-human experiments: grep fragments corresponding to species analyzed
## needs to be done only once ; first check if completed
if [ ! -f $PATH_T7"/hg19_fragments.tsv" ]
then
cd $PATH_T7
echo $PATH_T7
gzip -cd fragments.tsv.gz | grep hg19_ > hg19_fragments.tsv
#head hg19_fragments.tsv

## add sample identifier into fragment data
awk 'BEGIN { OFS = "\t" } $14 = "Wee1i." $4' hg19_fragments.tsv > temp.tsv
#head temp.tsv

# reorganize columns and remove hg19_ prefix from chr names
awk -v OFS="\t" '{print $1, $2, $3, $6, $5}' temp.tsv | sed 's/hg19_chr/chr/' > hg19_fragments.tsv
head hg19_fragments.tsv
rm temp.tsv
fi

## -----------------------

# same for ctrl
if [ ! -f $PATH_T8"/hg19_fragments.tsv" ]
then
cd $PATH_T8
echo $PATH_T8
gzip -cd fragments.tsv.gz | grep hg19_ > hg19_fragments.tsv
#head hg19_fragments.tsv

awk 'BEGIN { OFS = "\t" } $14 = "ctrl." $4' hg19_fragments.tsv > temp.tsv
#head temp.tsv

awk -v OFS="\t" '{print $1, $2, $3, $6, $5}' temp.tsv | sed 's/hg19_chr/chr/' > hg19_fragments.tsv
head hg19_fragments.tsv
rm temp.tsv
fi

if [ -f $PATH_T8"/hg19_fragments.tsv" ]
then
echo "fragment files with only fragments from species analyzed have been generated"
fi

## -----------------------
## merge fragments across treatments
## https://satijalab.org/signac/articles/faq.html?q=merge#how-do-i-merge-objects-with-signac
if [ ! -f $WD"/T7.T8.fragments.tsv.gz" ]
then
echo "merge fragments across treatments"
# merge files (avoids having to re-sort)
sort -m -k 1,1V -k2,2n $PATH_T7"/hg19_fragments.tsv" $PATH_T8"/hg19_fragments.tsv" > $WD"/T7.T8.fragments.tsv"

# block gzip compress the merged file
# -@ 4 uses 4 threads
bgzip -@ 4 $WD"/T7.T8.fragments.tsv"

# index the bgzipped file
tabix --preset=bed $WD"/T7.T8.fragments.tsv.gz"
fi

if [ -f $WD"/T7.T8.fragments.tsv.gz" ]
then
echo "merged fragment file has been generated"
fi

## -----------------------

## process with Signac and Seurat

## preprocess samples (extract human cells, create Seurat objects)
R CMD BATCH 01_Preprocessing.R

## merge peaks from both samples so that can create a peak count matrix for the merged object
cat temp_ctrl_peaks.txt | awk -v OFS='\t' -F':' '{print $1 OFS $2}' |awk -v OFS='\t' -F'-' '{print $1 OFS $2}' > temp_ctrl_peaks.bed

cat temp_Wee1i_peaks.txt | awk -v OFS='\t' -F':' '{print $1 OFS $2}' |awk -v OFS='\t' -F'-' '{print $1 OFS $2}' > temp_Wee1i_peaks.bed

cat temp_ctrl_peaks.bed temp_Wee1i_peaks.bed | sort -k1,1 -k2,2n > temp.sorted.bed
bedtools merge -i temp.sorted.bed > merged_peaks.bed

## keep peaks with min length 10
awk '{len = $3 - $2; if (len > 10) { print $_}}' merged_peaks.bed > merged_peaks_sizeFilt.bed

## region sizes before
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1     385     662    1079    1325   82503

## region sizes after merge
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1     371     656    1090    1318   82535

## calculate QC metrics and leave out poor quality cells
R CMD BATCH 02_qc-metrics.R

## normalize, cluster, find neighbors and plot the samples separately in lower dimension (UMAP)
R CMD BATCH 03_norm_dimred.R

## calculate gene activity matrix and save an object with this as default assay (ACTIVITY)
## the UMAP here is re-calculated using top genes; use this object to integrate with RNA-seq
R CMD BATCH 04_geneactivitymatrix.R

## calculate custom enhancer activity matrix and save an object with this as default assay
## the UMAP here is re-calculated using top enhancers
R CMD BATCH 05_enhanceractivitymatrix.v2.R
