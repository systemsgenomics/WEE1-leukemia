### Cell isolation and sample loading
RS4;11 and Nalm6

Dead Cell Removal Kit (#130-090-101, MACS miltenyi Biotech) was used to remove dead cells
- pools of human and mouse (mouse data not part of this study) cell lines were loaded together

### Library preparation 
ChromiumTM Chromium Single Cell 3´Reagent Kits v3 User guide CG000184 Rev A and libraries constructed using the 10X Genomics Chromium technology.
- Loading concentrations were 1700 cell/ul, 1000cell/ul for RS4;11 DMSO, AZD1775 cells; 500 cell/ul, 1300 cell/ul Nalm-6 DMSO, AZD1775, respectively.
- 
### Sequencing 
NovaSeq S2, aiming 50 000 reads per cell

### Preprocessing
cellranger_rna_script.sh

### Downstream analysis

human cell line data, 24 h
- SEURAT_singleCell_RNAseq_RS411_24h.R (QC, filtering, normalization; produces clustering used for rna-seq cell state assignment)
- scRNAseq_RS411_clusterDE_heatmap.R (marker genes for cell states)
- RS411_run_velocyto_score_genesets.py (produces UMAPs separated by treatment shown in figs)

- SEURAT_singleCell_RNAseq_Nalm6_24h.R (QC, filtering, normalization)
- integration_with_scRNAseq_24h_RS411_Nalm6.R (label transfer)
- Nalm6_run_velocyto_score_genesets.py (produces UMAPs separated by treatment shown in figs)

human cell line data, 10 h
- integration_with_scRNAseq_RS411.24h_RS411.10h_DMSOcellsonly.R (label transfer)
- integration_with_scRNAseq_RS411.24h_RS411.10h_Wee1icellsonly.R (label transfer)

human primary cell data
- SEURAT_singleCell_RNAseq_MLL7treated.R (QC, filtering, normalization)
- viz_mouse_sc_data_MLL7.R (marker gene plots)

- SEURAT_singleCell_RNAseq_MEF2D.R (QC, filtering, normalization)
- viz_10X_data_MEF2D.R (marker gene plots)

