### Scripts to preprocess GRO-seq data

The data from each sequencing experiment was run in two steps
1. prep
- main script: bixPREP.sh
- parameters are specified in prep_runconfig.txt
- files are specified in pdata.txt
- before run create folder prepres and subfolders trimconfigs and QC
- under the folder with scripts create a subfolder logs
- after run, check QC plots and specify possible filter/trim settings in trimconfig file created

2. gro
- main script: bixGRO-toTagdirs.sh
- parameters are specified in gro_testUser.txt
- before run retrieve mappability tracks and blacklisted regions: create_mappability.sh
- before run retrieve reference genome folder structure from iGenomes
- before run create folder gropres and subfolders aligned, built_indexes, tagdirs, intermed_files and QC

