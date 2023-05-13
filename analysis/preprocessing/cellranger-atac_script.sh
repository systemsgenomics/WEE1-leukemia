#!/bin/bash

cellranger_atac=cellranger-atac-1.2.0/cellranger-atac
ref=refdata-cellranger-atac-hg19-and-mm10-1.2.0
fastq_path=data/scATACseq/raw_data/fastq/
mem=360
cores=40

for i in {1..8}
do
        id=T$i

        $cellranger_atac count --id=$id --reference=$ref --fastqs=$fastq_path$id --sample=$id --localcores=$cores --localmem=$mem --disable-ui
done
