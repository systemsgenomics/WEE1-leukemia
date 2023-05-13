#!/bin/bash

cellranger_multi=cellranger-arc-2.0.1/cellranger-arc
ref=refdata-cellranger-arc-GRCh38-2020-A
libraries=data/scMultiome/libraries.csv
mem=200
cores=24



$cellranger_multi count --id=$id --reference=$ref --libraries=$libraries --localcores=$cores --localmem=$mem --disable-ui
