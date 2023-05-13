## sh $SCRIPTS$(echo homer/tagdirs.sh $PROJ $WWW $REF $PUBLIC $FLENGTH $CUSTOM $GVER)

PROJ=$1
WWW=$2
REF=$3
PUBLIC=$4
FLENGTH=$5
CUSTOM=$6
GVER=$7
URL=$8
PIPE=$9

DESTDIR=$PIPE"/grores/"


replace=""

cd $DESTDIR"aligned/"$PROJ"/"
echo " "
echo "creating Homer tagdirs"
date
echo " "
echo "processing these alignment files:"
ls | grep ".bed"
for i in $( ls | grep ".bed"); do
	name=$( echo $i | sed -e "s/.bed/$replace/g")
	echo "SAMPLE NAME: "$name
	
	echo "removing reads that overlap blacklisted regions (encode no strand info, custom stranded) or snoRNA loci"
	## exclude reads that overlap blacklisted regions, naming will be the same whether or not this is done or not
	if [ -f $CUSTOM"mappable/nonmappable_coords/"$REF"/"$REF"-blacklist.bed" ]
	then
	echo $i
	echo intersectBed -a $i -b $CUSTOM"mappable/nonmappable_coords/"$REF"/"$REF"-blacklist.bed" -f 0.99 -v > $name"_filtered.bed"
	intersectBed -a $i -b $CUSTOM"mappable/nonmappable_coords/"$REF"/"$REF"-blacklist.bed" -f 0.99 -v > $name"_filtered.bed"

		if [ -f $CUSTOM$REF"del.bed" ]
		then
		mv $name"_filtered.bed" "tempX.bed"
		## remove reads mapping to custom blacklisted or snoRNA loci strand-specific 
		intersectBed -a "tempX.bed" -b $CUSTOM$REF"del.bed" -f 0.99 -v -s > $name"_filtered.bed"
		rm "tempX.bed"
		else
		echo "could not perform custom blacklist filtering, file not found"
		echo $CUSTOM$REF"del.bed"
		fi

	else
	echo "could not perform encode blacklist filtering, file not found"
	echo $CUSTOM"mappable/nonmappable_coords/"$REF"/"$REF"-blacklist.bed"
	mv $i $name"_filtered.bed"
	fi



	makeTagDirectory $DESTDIR"tagdirs/"$name $name"_filtered.bed" -genome $REF -checkGC -fragLength $FLENGTH
	makeMultiWigHub.pl $name $REF -url $URL -webdir $PIPE"/hubs/"$REF"/groseq/" -strand -d $DESTDIR"tagdirs/"$name -force
	cp -r $PIPE"/hubs/"$REF"/groseq/"$name $WWW
	chmod -R 775 $WWW"/"$name
	cd $DESTDIR"tagdirs/"$name
	R CMD BATCH $PIPE"/gro/scripts/qc/plotQC.R"
	if [ $PUBLIC = "FALSE" ]
	then
	cp -r $DESTDIR"tagdirs/"$name $PIPE"/toIDA/tagdirs"
	fi


	cd $DESTDIR"aligned/"$PROJ"/"
	

	date
done



