cpus=8
mem=25000
gtf=refdata-cellranger-hg19-and-mm10-3.0.0/genes/genes.gtf
cellranger_output= # provide path here

velocyto run10x -v -@ $cpus --samtools-memory $mem $cellranger_output $gtf
