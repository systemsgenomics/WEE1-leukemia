## Solexa preprocessing pipeline
## Required software installations:
## R (standard libraries)
## Homer
## FastQC
## NGSQCToolkit
## (assumes that these software are found on the path)

## USAGE: nohup sh ./bixPREP.sh /path/to/runconfigs/myparam_file.txt &
## for example nohup sh ./bixPREP.sh $HOME/prep/runconfigs/preprocess_template.txt &

## what have you done before:
## 1. create pdata.txt
## 2. saved to where your fastq/sra files are located
## 3. created run configuration file (saved under ../runconfigs)

## The purpose of pdata.txt and parameter configuration file (see below) is to save YOU time to create dirs and move files or set parameters 
## ->  you just need to describe them at start
## TIP: if seq was done in multiple runs, the pipeline will later combine them to one file
## if you set identifical gsmid in pdata.txt !

## You will know how far the script is running by checking the logs folder and the run specific log file

## OUTPUT:
## re-named files
## folder to be copied to backups (IDA)
## QC plots
## trimconfig files that can be used to
## 1) further trim reads based on QC plot, optional step before next pipeline
## 2) specify whether to include all runs of the same sample to analysis (useful if e.g. one sequencing failed to omit it)

##----------------------------------------------------------
## Destdir specifies the pipeline output folder structure path where all your processed files will end up 
## All results will be going to folders named using project_name based on pdata entries project, assay and GSEid


##**********************************************************
## PARAMETERS THAT ARE SET BASED ON RUN CONFIG

PARAMS=$1
. $PARAMS

echo "parameters"
echo "data is public: " $PUBLIC
echo "name of the user running analysis: " $USER
echo "folder name for files used by pipeline run: " $FOLDER
echo ".sra/fastq input file location: " $FILES
echo "pipeline folder: " $PIPE
echo "data transferred after processing to: " $PIPE"/prepres/"
echo "files to be kept also copied under: " $PIPE"/toIDA/"$USER 

echo "file type [sra/fastq] : " $TYPE
echo "data needs splitting by barcodes : " $PROCESS_BARCODES
echo "quality score offset :" $OFFSET
echo "sample is :" $SEQWAY"-end"
echo "polyA trimming needed :" $TRIM
echo "polyA strech to use in trimming: " $TRIMA
echo "read length cutoff :" $TRIMLENGTH

## the pipeline will consider this one project, named as specified by FOLDER
PROJ=$FOLDER
DESTDIR=$PIPE"/prepres/"

##**********************************************************

## PRE-SET PARAMETERS (TYPICALLY NO NEED TO CHANGE)

## Main directory of the pipeline. This should be located in your home folder. Contains scripts folder. cd to that folder to run this script
PIPELINE=$(echo /home/groups/biowhat/ngs_pipelines/prep/)

## tool lists, these are on common path

## http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=std
##SRATOOLS

## http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/
##FASTQC

## http://biowhat.ucsd.edu/homer/download.html
##HOMER

# other paths
LOGS=$PIPELINE$(echo logs/)
SCRIPTS=$PIPELINE$(echo scripts/)

##**********************************************************

## PARAMETERS DETERMINED DURING RUN

## the process ID and time stamp are used to create logs and run-specific naming
PID=$(echo $$)
TIME=$(date +"%F_%H%M")


##**********************************************************
##**********************************************************
##----------------------------------------------------------
##  run info

echo "STARTING ANALYSIS IN FOLDER" $FILES
echo "YOUR PROCESS ID IS "
echo $PID

echo "A LOG OF ALL SCRIPTS EXECUTED WILL BE STORED AT"
echo $LOGS$PROJ"_"$TIME$PID.log
echo "ANALYSIS LOG  $LOGS$PROJ"_"$TIME$PID.log" > $LOGS$PROJ"_"$TIME$PID.log

date 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
echo . $PARAMS  2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log

## sanity check
if [ -d $DESTDIR ]
then
echo "project folder will appear in pre-existing" $DESTDIR 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
else
mkdir $DESTDIR
echo "project folder will appear in directory created:" $DESTDIR 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
fi

## clean-up for fresh start
if [ -d $DESTDIR$PROJ ]
then
rm -r $DESTDIR$PROJ
echo "removed pre-existing" $DESTDIR$PROJ 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
fi

if [ -d $PIPE"/toIDA/"$USER"/"$FOLDER ]
then
rm -r $PIPE"/toIDA/"$USER"/"$FOLDER
echo "removed pre-existing" $PIPE"/toIDA/"$USER"/"$FOLDER 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
fi

## make folders to keep run files and IDA backup ready files
mkdir $DESTDIR$PROJ
echo "converted fastq files of this project will appear in directory created:" $DESTDIR$PROJ 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log

mkdir $PIPE"/toIDA/"$USER
mkdir $PIPE"/toIDA/"$USER"/"$FOLDER
echo "files for IDA backups will appear in directory created:" $PIPE"/toIDA/"$USER"/"$FOLDER 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
##----------------------------------------------------------

## file re-naming based on pdata and moving into pipeline folder structure (prepres)
echo "uncompressing and sorting files" 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
cd $FILES

## make a copy for IDA transfer if this is unpublished data
if [ "$PUBLIC" = "FALSE" ]
then
## files should be compressed
cp *.gz $PIPE"/toIDA/"$USER"/"$FOLDER
## in case not...will also try to put there fastq and txt files (not bam files that sometimes come from seq core)
cp *.fastq $PIPE"/toIDA/"$USER"/"$FOLDER
cp *.txt $PIPE"/toIDA/"$USER"/"$FOLDER
fi

## if there is anything compressed, get it ready for run by uncompressing
gunzip *.gz

if [ "$SEQWAY" = "paired" ]
then
R --slave --vanilla --file=$SCRIPTS$(echo preparePE.R) --args $OFFSET $TYPE $DESTDIR $PROCESS_BARCODES 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
else
	if [ "$PROCESS_BARCODES" = "FALSE" ]
	then
	R --slave --vanilla --file=$SCRIPTS$(echo prepare.R) --args $OFFSET $TYPE $DESTDIR $TRIM $TRIMA 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
	else
	R --slave --vanilla --file=$SCRIPTS$(echo prepare_split_by_barcode.R) --args $OFFSET $TYPE $DESTDIR $TRIM $TRIMA 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
	fi
fi
date 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log



## -------------- TRIMMING and QC ------------------------------------------------------------------
## Single-end reads:
## 1) Trim polyA -> remove too short reads (based on TRIMLENGTH) -> perform QC analysis
## 2) Just remove too short reads (useful if reads differ in length, not too useful otherwise)
## 3) QC reads after filtering, or if no filtering, then QC anyway 
## 4) copy reads that are now ready to be submitted to GEO to the ngs_pipelines/toIDA folder

## Paired-end reads:
## previous step already produced QC reports and made copy under ngs_pipelines/toIDA folder
## trimming based on polyA or length performed during the actual seq run
 
echo " "
echo "TRIMMING" 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
date
echo " "

cd $DESTDIR$PROJ
replace=""


if [  "$SEQWAY" = "single" ]
then
	 if [ "$TRIM" = "TRUE" ]
	 then
	 echo "starting trimming"
		 ## just trim based on length, a precautionary step which you should do only if needed
		 if [ "$TRIMA" = "NA" ]
		 then
		    echo "just length trimming"
		    homerTools trim -min $TRIMLENGTH *.fastq* 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
		    date 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
		    replace=""
		    for i in $( ls | grep ".trimmed"); do
			     $FASTQC$(echo fastqc) $i --outdir=$DESTDIR"QC/"  2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
			    if [ "$PUBLIC" = "FALSE" ]
			    then			
			    cp $i $PIPE"/toIDA/"$USER"/"$FOLDER
			    gzip $PIPE"/toIDA/"$USER"/"$FOLDER/$i
			    fi			
			    cp $DESTDIR"QC/"$i"_fastqc.zip" $PIPE"/toIDA/"$USER"/"$FOLDER
			    date 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
		    done
		 else
		    echo "polyA and length trimming"
		    ## typical trimming based on polyA, followed by length filter
		    homerTools trim -3 $TRIMA  *.fastq* 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
		    ##will add ending *.trimmed
		    date 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
		    ## Trim here the reads to be at least X bp long, set by $TRIMLENGTH parameter
		    homerTools trim -min $TRIMLENGTH *.trimmed 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
		    date 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log	
		    ##will add ending *.trimmed
		    for i in $( ls | grep ".trimmed.trimmed"); do
			    $FASTQC$(echo fastqc) $i --outdir=$DESTDIR"QC/"  2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
			    if [ "$PUBLIC" = "FALSE" ]
			    then			
			    cp $i $PIPE"/toIDA/"$USER"/"$FOLDER
			    gzip $PIPE"/toIDA/"$USER"/"$FOLDER/$i
			    fi
			    cp $DESTDIR"QC/"$i"_fastqc.zip" $PIPE"/toIDA/"$USER"/"$FOLDER
			    date 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
		    done
		 fi
    fi

	if [ "$TRIM" = "FALSE" ]
	then
	echo "no trimming performed, check sample QC"
	## this will happen if no trimming done
		for i in $( ls | grep ".fastq"); do
		$FASTQC$(echo fastqc) $i --outdir=$DESTDIR"QC/"  2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
			if [ "$PUBLIC" = "FALSE" ]
			then
			cp $i $PIPE"/toIDA/"$USER"/"$FOLDER
			gzip $PIPE"/toIDA/"$USER"/"$FOLDER/$i
			fi
			replace=""
			new=$( echo $i | sed -e "s/.fastq/$replace/g")
			cp $DESTDIR"QC/"$new"_fastqc.zip" $PIPE"/toIDA/"$USER"/"$FOLDER
		done	
	fi

fi

date 2>&1 | tee -a $LOGS$PROJ"_"$TIME$PID.log
## QC done 
## edit the respective *trimconfig.txt files to set trimming options
##----------------------------------------------------------

echo "copy run info under ngs_pipelines/toIDA/userName/folderName"
cp $PARAMS $PIPE"/toIDA/"$USER"/"$FOLDER 
cp $FILES"/pdata.txt" $PIPE"/toIDA/"$USER"/"$FOLDER
cp $LOGS$PROJ"_"$TIME$PID.log $PIPE"/toIDA/"$USER"/"$FOLDER

echo "prep done"


