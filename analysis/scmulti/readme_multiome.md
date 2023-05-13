### Nuclei isolation and sample loading

10x Genomics protocol (User guide CG000365 RevB), optimized for the RS4;11 cell line with a 3-minute cell lysis time.
- Loading concentrations 7200 nuclei/ul for DMSO cells and 4860 nuclei/ul for AZD1775 cells

### Library preparation

Chromium Next GEM Single Cell Multiome ATAC + Gene Expression User guide CG000338 Rev E

### Sequencing

ATAC: NovaSeq PE50 aiming 20 000 reads per cell depth
RNA: NovaSeq PE150 aiming 50 000 reads per cell depth

### Preprocessing

cellranger-multiome_script.sh

### Downstream analysis

- scATAC data processed to generate similar outputs that can be compared to 24 h samples

main script: main_atac.sh (extract human fragment data)

Analysis with Signac and Seurat R packages

Seurat: Stuart, Butler et al. 2019, version 3.1.1 and Signac: Stuart, Butler et al. 2019, version 0.1.6

run separately in this order:

01_Preprocessing.R
- writes out bed files of peaks from DMSO and AZD1775 samples bedtools merge (part of main script)

01.2_Preprocessing.R

02_qc-metrics.R

03_norm_dimred.R

