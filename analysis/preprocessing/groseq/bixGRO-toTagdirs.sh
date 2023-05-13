## GRO-seq (Homer-based) pipeline
## assumes that data has been preprocessed using the prep pipeline

## USAGE: sh ./bixGRO.sh /path/to/myparam_file.txt



##----------------------------------------------------------

## The parameter configuration takes care of setting paths to where your data is located, options for file conversion and alignment
## Set the pipeline folder as destdir where all your processed files will end up

## PARAMETERS THAT ARE SET BASED ON RUN CONFIG myparam_file.txt


PARAMS=$1
. $PARAMS

echo "parameters"
echo "data is public: " $PUBLIC
echo "name of the user running analysis: " $USER
echo "folder name for files used by pipeline run: " $FOLDER
echo "analysis step : " $STEP
echo "perform just this step : " $JUSTONE
echo "pipeline folder: " $PIPE
echo "will look for preprocessed fastq input files in: " $PIPE"/prepres/"
echo "data transferred after processing to: " $PIPE"/grores/"
echo "quality score offset :" $OFFSET
echo "% nt in read passing quality cutoff : " $QP
echo "quality score cutoff :" $QMIN
echo "genome from :" $GVER
echo "reference genome :" $REF
echo "max mismatches allowed (bowtie v) :" $BOWTIEv
echo "max number of hits in genome (bowtie m) :" $BOWTIEm
echo "report this many hits if many found (bowtie k) :" $BOWTIEk
echo "fragment length estimate given :" $FLENGTH
echo "stack size is defined to be :" $PC
echo "number of cores to parallelize the run to: " $PARALLEL
echo "clean-up intermediate files: " $CLEAN
echo "URL for hubs: " $URL

PROJ=$FOLDER
FILES=$PIPE"/prepres/"
DESTDIR=$PIPE"/grores/"


## Main pipeline folder  

PIPELINE=$PIPE
# other paths
LOGS=$PIPELINE$(echo logs/)
SCRIPTS=$PIPELINE
CUSTOM=$PIPE"/custom_coords/pipefiles/"
WWW=$(echo /home/www/hubs/)$REF"/groseq/"
echo "path for hubs: " $WWW

## PARAMETERS THAT ARE PRE-SET (no need to modify typically)
## Genome file locations

GPATH=$(echo /home/work/public/iGenomes/)

if [ $REF = "hg19" -o $REF = "hg18" ]
then
GSPEC=$(echo Homo_sapiens)
fi

if [ $REF = "mm9" -o $REF = "mm10" -o $REF = "mm8" ]
then
GSPEC=$(echo Mus_musculus)
fi

if [ $REF = "rn4" ]
then
GSPEC=$(echo Rattus_norvegicus)
fi

echo "if species is not human, mouse or rat (hg18, hg19, mm8, mm9, mm10, rn4) the run will fail"

## required index files for bowtie and samtools
BgenINDEX=$(echo $GPATH$GSPEC"/"$GVER"/"$REF"/Sequence/BowtieIndex/")
BdelINDEX=$(echo $GPATH$GSPEC"/"$GVER"/"$REF"/Sequence/AbundantSequences/")


SINDEX=$(echo $GPATH$GSPEC"/"$GVER"/"$REF"/Sequence/WholeGenomeFasta/")
ANNODIR=$(echo $GPATH$GSPEC"/"$GVER"/"$REF"/Annotation/Genes/")


## tools, most of these are on common path, then no need to call with full path

## http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=std
#SRATOOLS

## http://hannonlab.cshl.edu/fastx_toolkit/
##FASTX

## http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/
##FASTQC

##http://bowtie-bio.sourceforge.net/index.shtml
##BOWTIE

## if Homer is not on your path, you can add this line to point to homer directory with binaries
##HOMER

## UCSC binaries, deprecated
#UCSCbin=$(echo /home/bgrp3/merjahe/ucsc_binaries/)
## SAM tools
## http://samtools.sourceforge.net/
##SAMTOOLS=$(echo "")
## requires also index of genome: it is found under iGenomes annotations, no need to do anything but if you want your own see below



##**********************************************************

## PARAMETERS DETERMINED DURING RUN

## the process ID and time stamp are used to create logs and run-specific naming
PID=$(echo $$)
TIME=$(date +"%F_%H%M")

##----------------------------------------------------------


##**********************************************************
##**********************************************************
##----------------------------------------------------------
## THESE STEPS TAKE YOU FROM PREPROCESSED READS TO ALIGNED READS
echo "STARTING ANALYSIS IN FOLDER" $FILES
echo "YOUR PROCESS ID IS "
echo $PID

echo "A LOG OF ALL SCRIPTS EXECUTED WILL BE STORED AT"
echo $LOGS$PROJ"_"$SAMPLE"_"$TIME$PID.log
echo "ANALYSIS LOG  $LOGS$PROJ"_"$TIME$PID.log" > $LOGS$PROJ"_"$TIME$PID.log

echo . $PARAMS 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log

if [ -f $DESTDIR"built_indexes/"$GVER$REF"_seq2del.1.ebwt" ]
then
echo "abundant sequence index found" 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
else
echo " will build index of abundant (contamination) sequences to your destdir/built_indexes directory" 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
cd $DESTDIR
mkdir built_indexes
cd built_indexes
cp $BdelINDEX*.fa .
## add any custom seqs to remove to custom_coords, will be added to index
## these should be seqs that relate to library prep (not species-specific)
cp $CUSTOM*.fa .
## not using polyC removal
rm polyC.fa
replace=""
echo "building index from" 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
ls -m *.fa 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
fastas=$(ls -m *.fa)
toindex=$( echo $fastas | sed -e "s/ /$replace/g")
## create on the fly the abundant seq bowtie index
bowtie-build -f $toindex $GVER$REF"_seq2del" 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
rm *.fa
fi

if [ $STEP = "1" ]
then

	## -------------- (OPTIONAL) TRIMMING (FASTX) AND FILTERING-------------------------------------------------

	echo " "
	echo "FILTER BY QUALITY, PERFORM EXTRA TRIMMING" 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
	date 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
	echo " "

	R --slave --vanilla --file=$SCRIPTS$(echo qc/filter_trim.R) --args $FILES $PROJ $OFFSET $QP $QMIN $DESTDIR 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log

	if [ $CLEAN = "TRUE" ]
	then
	rm -r $FILES"/"$FOLDER
	fi

## are we done or should we proceed to next steps?
	if [ "$JUSTONE" = "TRUE" ]
	then
	STEP="0"
	else
	echo "filtering done" 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
	STEP="2"
	fi



fi



if [ $STEP = "2" ]
then
echo " "
	## -------------- GENOME ALIGNMENT (BOWTIE)-------------------------------------------------------------
	echo " "
	echo "START ALIGNMENTS" 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
	date 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
	echo " "

 
	sh $SCRIPTS$(echo bowtie/justBowtie_fq.sh $PROJ  $DESTDIR $PARALLEL $REF $BOWTIEv $BOWTIEm $BOWTIEk $BgenINDEX $SINDEX $GVER) 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log

	if [ $PUBLIC = "FALSE" ]
	then
	mkdir $PIPE"/toIDA/"$USER"/"$FOLDER"/alignments/"
	cp -r $DESTDIR"/aligned/"$FOLDER"/QC/" $PIPE"/toIDA/"$USER"/"$FOLDER"/alignments/"
	cp $DESTDIR"/aligned/"$FOLDER"/"*sorted* $PIPE"/toIDA/"$USER"/"$FOLDER"/alignments/"
	gzip $DESTDIR"/aligned/"$FOLDER"/"*.bowtie
	mv $DESTDIR"/aligned/"$FOLDER"/"*.bowtie.gz $PIPE"/toIDA/"$USER"/"$FOLDER"/alignments/"
	fi

	if [ $CLEAN = "TRUE" ]
	then
	rm -r $DESTDIR"/intermed_res/"$FOLDER
	fi


## are we done or should we proceed to next steps?
	if [ "$JUSTONE" = "TRUE" ]
	then
	STEP="0"
	else
	echo "alignments done" 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
	STEP="3"
	fi



fi




if [ $STEP = "3" ]
then

	## -------------- CUSTOM TAGDIRS -------------------------------------------------

	echo " "
	echo "CREATING CUSTOM TAGDIRS AND VIZ" 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
	date 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
	echo " "


	sh $SCRIPTS$(echo homer/tagdirs_custom.sh $PROJ $WWW $REF $PUBLIC $FLENGTH $CUSTOM $GVER $URL $PIPE $PUBLIC) 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log

	if [ $CLEAN = "TRUE" ]
	then
	rm -r $DESTDIR"/aligned/"$FOLDER
	fi	

	echo "tagdirs done" 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
fi


echo "copy run info under toIDA/userName/folderName"
cp $PARAMS $PIPE"/toIDA/"$USER"/"$FOLDER 
cp $FILES"/trimconfigs/"$FOLDER"trimconfig.txt" $PIPE"/toIDA/"$USER"/"$FOLDER
cp $LOGS$PROJ"_"$TIME$PID.log $PIPE"/toIDA/"$USER"/"$FOLDER

 
