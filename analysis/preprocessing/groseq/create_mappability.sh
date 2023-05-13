## Create mappability bed files from UCSC table browser mappability tracks
## Will also include ENCODE blacklisted regions as non-mappable (reads filtered)

##sh $SCRIPTS$(echo custom/create_mappability.sh $REF $CUSTOM)

# genome version
REF=$1
# output folder
CUSTOM=$2
# download yes: 1, no: 0
DL=$3

mkdir $CUSTOM"mappable"
mkdir $CUSTOM"mappable/nonmappable_coords"

## go to this folder first
cd $CUSTOM"mappable/nonmappable_coords"
mkdir $REF
cd $REF

## files downloaded from UCSC, mappability tracks, wget downloads them all
if [ $DL = "1" ]
then
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/$REF/encodeDCC/wgEncodeMapability/*
##wget ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeMapability/*

	for i in $(echo 36 40 50 75 100) ;do
		echo bigWigToBedGraph "wgEncodeCrgMapabilityAlign"$i"mer.bigWig" $REF"_wgEncodeCrgMapabilityAlign"$i"mer.bedGraph"
		bigWigToBedGraph "wgEncodeCrgMapabilityAlign"$i"mer.bigWig" $REF"_wgEncodeCrgMapabilityAlign"$i"mer.bedGraph"
	done
gzip *.bigWig


## then retrieve blacklisted regions documented here https://sites.google.com/site/anshulkundaje/projects/blacklists

	if [ $REF = "mm9" ]
	then
	wget http://www.broadinstitute.org/~anshul/projects/mouse/blacklist/mm9-blacklist.bed.gz
	gunzip mm9-blacklist.bed.gz
	mv mm9-blacklist.bed xx.bed
	awk '{print $0"\tencode_blacklisted\t.\t."}' xx.bed > $REF"-blacklist.bed"
	rm xx.bed
	fi


	if [ $REF = "hg19" ]
	then
	## files downloaded include also blacklisted regions 
	gunzip wgEncodeDacMapabilityConsensusExcludable.bed.gz
	mv wgEncodeDacMapabilityConsensusExcludable.bed $REF"-blacklist.bed"
	fi

fi

##  all kmers

for i in $(echo 36 40 50 75 100) ;do

	gunzip $REF"_wgEncodeCrgMapabilityAlign"$i"mer.bedGraph.gz"
	# Find those that do not match uniquely enough, has over 2 matches
	awk '$4<=0.33 {print $0}' $REF"_wgEncodeCrgMapabilityAlign"$i"mer.bedGraph" | cut -f 1-3 | awk '{print $0"\t.\t.\t."}' > $REF"_wgEncodeCrgMapabilityAlign"$i"mer_nonmappable_temp.bed"
	gzip $REF"_wgEncodeCrgMapabilityAlign"$i"mer.bedGraph"
	cat $REF"-blacklist.bed" $REF"_wgEncodeCrgMapabilityAlign"$i"mer_nonmappable_temp.bed" | cut -f 1-5 > $REF"_wgEncodeCrgMapabilityAlign"$i"mer_nonmappable_temp2.bed"

	rm $REF"_wgEncodeCrgMapabilityAlign"$i"mer_nonmappable_temp.bed"

	## need to sort the file
	sort -k1,1 -k2,2n $REF"_wgEncodeCrgMapabilityAlign"$i"mer_nonmappable_temp2.bed" > $REF"_wgEncodeCrgMapabilityAlign"$i"mer_nonmappable_temp2.bed.sorted.bed"

	rm $REF"_wgEncodeCrgMapabilityAlign"$i"mer_nonmappable_temp2.bed"

	##and then merge, adds extra tab at the end so that it is (R-)easy to just add strand info if needed
	bedtools merge -i $REF"_wgEncodeCrgMapabilityAlign"$i"mer_nonmappable_temp2.bed.sorted.bed" > "../"$REF"_wgEncodeCrgMapabilityAlign"$i"mer_nonmappable.bed"

	rm $REF"_wgEncodeCrgMapabilityAlign"$i"mer_nonmappable_temp2.bed.sorted.bed"
done
cd ..

