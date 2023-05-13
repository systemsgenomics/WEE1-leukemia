
##sh $SCRIPTS$(echo bowtie/justBowtie_fq.sh $PROJ  $DESTDIR $PARALLEL $REF $BOWTIEv $BOWTIEm $BOWTIEk $BgenINDEX $SINDEX $GVER)

PROJ=$1
DESTDIR=$2
PARALLEL=$3
REF=$4
BOWTIEv=$5
BOWTIEm=$6
BOWTIEk=$7
BgenINDEX=$8
SAMTOOLS=${9}
GVER=${10}

## -------------- GENOME ALIGNMENT (BOWTIE)-------------------------------------------------------------


echo " "
echo "START ALIGNMENTS FOR FILTERED READS"
date
echo " "

cd $DESTDIR"intermed_res/"$PROJ
#echo "removing existing intermediate results, if none exist don't be alarmed by No such file or directory msg that follows"
rm cleaned*
rm *.bowtie

mkdir rRNA_res

ls | grep "qfilt_"
## alignment to rRNA, adapter sequence, polyA, polyC, etc..
for i in $( ls | grep "qfilt_"); do
	echo "removing abundant (uninteresting) sequences from : "$i
	bowtie $DESTDIR"built_indexes/"$GVER$REF"_seq2del" -q $i $i"_hs_rRNA".bowtie --best --strata -v 1 -m 3 -k 1 -p $PARALLEL --un "cleaned_"$i
	date
done


if [ -d $DESTDIR"aligned/"$PROJ ]
then
echo "removing existing alignments for this project"
rm -r $DESTDIR"aligned/"$PROJ
fi

mkdir $DESTDIR"aligned/"$PROJ
mkdir $DESTDIR"aligned/"$PROJ"/QC"
## alignment of fastq files, outputs bowtie format file (compatible with Homer)
ls | grep "cleaned"
for i in $( ls | grep "cleaned"); do
	echo "genome alignment and bam conversion"
	echo "FILENAME: "$i
	clean1=$( echo $i | sed -e "s/.fastq/$replace/g")
	cleanend=$( echo $clean1 | sed -e "s/.trimmed/$replace/g")
	cleanstart=$( echo $cleanend | sed -e "s/rRNAcleaned_qfilt_/$replace/g")
	name=$( echo $cleanstart | sed -e "s/cleaned_qfilt_/$replace/g")
	echo "OUTPUT FILE: " $name"_"$GVER"_"$REF.bowtie
	echo bowtie $BgenINDEX"genome" -q $i $DESTDIR"aligned/"$PROJ"/"$name"_"$GVER"_"$REF.bowtie --best --strata -v $BOWTIEv -m $BOWTIEm -k $BOWTIEk -p $PARALLEL
	bowtie $BgenINDEX"genome" -q $i $DESTDIR"aligned/"$PROJ"/"$name"_"$GVER"_"$REF.bowtie --best --strata -v $BOWTIEv -m $BOWTIEm -k $BOWTIEk -p $PARALLEL
	echo "OUTPUT FILE: " $name"_"$GVER"_"$REF.sam
	bowtie $BgenINDEX"genome" -q $i $DESTDIR"aligned/"$PROJ"/"$name"_"$GVER"_"$REF.sam --best --strata -v $BOWTIEv -m $BOWTIEm -k $BOWTIEk -p $PARALLEL --sam
	date

## CONVERTING TO BED FOR FILTERING WITH BED TOOLS THEN PROCEEDING TO HOMER
	## convert to bed
	echo "converting to bed"
	awk -F"\t" '{print length($6)"\t"$0}' $DESTDIR"aligned/"$PROJ"/"$name"_"$GVER"_"$REF.bowtie | awk -F"\t" '{print $5+$1"\t"$0}' | awk -F"\t" '{print $5"\t"$6"\t"$1"\t"$3"\t"$9"\t"$4}' > $DESTDIR"aligned/"$PROJ"/"$name"_"$GVER"_"$REF.bed

	date


	
	## make also sam and convert to .bam for saving space (sam file deleted)
	## NOTICE THAT THESE FILES KEEP ALL READS, USE FLAGS TO RESOLVE USABLE READS
	samtools view -b -t $SINDEX"genome".fa.fai -o $DESTDIR"aligned/"$PROJ"/"$name"_"$GVER"_"$REF.bam $DESTDIR"aligned/"$PROJ"/"$name"_"$GVER"_"$REF.sam
	samtools sort $DESTDIR"aligned/"$PROJ"/"$name"_"$GVER"_"$REF.bam $DESTDIR"aligned/"$PROJ"/"$name"_"$GVER"_"$REF.sorted
	samtools index $DESTDIR"aligned/"$PROJ"/"$name"_"$GVER"_"$REF.sorted.bam

	samstat $DESTDIR"aligned/"$PROJ"/"$name"_"$GVER"_"$REF.bam
	mv $DESTDIR"aligned/"$PROJ"/"$name"_"$GVER"_"$REF.bam.samstat.html $DESTDIR"aligned/"$PROJ"/QC"
	rm $DESTDIR"aligned/"$PROJ"/"$name"_"$GVER"_"$REF.sam
	rm $DESTDIR"aligned/"$PROJ"/"$name"_"$GVER"_"$REF.bam
	date
done






