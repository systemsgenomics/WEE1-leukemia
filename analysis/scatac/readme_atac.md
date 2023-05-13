### Nuclei isolation and sample loading

RS4;11, 24 h timepoint
10x Genomics Nuclei Isolation for Single Cell ATAC Sequencing user guide CG000169 Rev D.
- pools of human and mouse (mouse data not part of this study) cell lines were loaded together

### Library preparation
Chromium Single Cell ATAC Library & Gel Bead and Chip E Single Cell ATAC Kit (10x Genomics), user guide CG000168 Rev C

### Sequencing  
NovaSeq S2, aiming 50 000 reads per nucleus

### Preprocessing
cellranger-atac_script.sh
- version 1.2.0

### Downstream analysis
main script: main_atac.sh
- extract human fragment data

### Analysis with Signac and Seurat R packages
Seurat: Stuart, Butler et al. 2019, version 3.1.1 and Signac: Stuart, Butler et al. 2019, version 0.1.6

included to main script or run separately in this order:

01_Preprocessing.R
- writes out bed files of peaks from DMSO and AZD1775 samples
- bedtools merge (part of main script)

01.2_Preprocessing.R

02_qc-metrics.R

03_norm_dimred.R

04_geneactivitymatrix.R

05_enhanceractivitymatrix.v2.R
