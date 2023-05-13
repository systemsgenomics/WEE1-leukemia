#!/bin/bash

cellranger_rna=cellranger-3.0.2/cellranger-cs
ref=refdata-cellranger-hg19-and-mm10-3.0.0
fastq_path=data/scRNAseq/raw_data/fastq/
mem=360
cores=40

for i in {1..8}
do
        id=S$i

        $cellranger_rna count --id=$id --transcriptome=$ref --fastqs=$fastq_path$id --sample=$id --localcores=$vores --localmem=$mem --disable-ui
done
