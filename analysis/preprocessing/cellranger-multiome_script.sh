#!/bin/bash

cellranger_multi=cellranger-arc-2.0.1/cellranger-arc
ref=refdata-cellranger-arc-GRCh38-2020-A
libraries=data/scMultiome/libraries.csv
mem=200
cores=24



$cellranger_multi count --id=$id --reference=$ref --libraries=$libraries --localcores=$cores --localmem=$mem --disable-ui

# only-ATAC-style, to compare to 24 h result
cellranger-atac-1.2.0/cellranger-atac-cs/1.2.0/bin/count --id=MA1-SI-NA-A3_HMLKTDRXY --reference=refdata-cellranger-atac-hg19-and-mm10-1.2.0 --fastqs=scATAC/MA1 --sample=MA1-SI-NA-A3_HMLKTDRXY --localcores=40 --localmem=360 --disable-ui
cellranger-atac-1.2.0/cellranger-atac-cs/1.2.0/bin/count --id=MA2-SI-NA-B3_HMLKTDRXY --reference=refdata-cellranger-atac-hg19-and-mm10-1.2.0 --fastqs=scATAC/MA2 --sample=MA2-SI-NA-B3_HMLKTDRXY --localcores=40 --localmem=360 --disable-ui

# only-RNA-style, to compare to 24 h result
cellranger-6.0.2/cellranger count --id=M1R-SCI7T014-SCI5T014_HK2L7DSX3 --transcriptome=refdata-cellranger-hg19-3.0.0 --fastqs=M1R --sample=M1R-SCI7T014-SCI5T014_HK2L7DSX3 --chemistry=ARC-v1 --localcores=24 --localmem=200
cellranger-6.0.2//cellranger count --id=M2R-SCI7T026-SCI5T026_HK2L7DSX3 --transcriptome=refdata-cellranger-hg19-3.0.0 --fastqs=M2R --sample=M2R-SCI7T026-SCI5T026_HK2L7DSX3 --chemistry=ARC-v1 --localcores=24 --localmem=200
