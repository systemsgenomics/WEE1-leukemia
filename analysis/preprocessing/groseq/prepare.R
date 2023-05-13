## version that considers the case that several files should be combined to form one sample

##R --slave --vanilla --file=$SCRIPTS$(echo prepare.R) --args  $OFFSET $TYPE $DESTDIR $TRIM $TRIMA
args=(commandArgs(TRUE))
#print(args)


offset=args[1]
type=args[2]
destdir=args[3]
trim=as.logical(args[4])
trima=as.logical(trim&(args[5]!="NA"))

cwd=getwd()
fileCombine=F

## information about files to be processed during the run (use template and take care to fill it carefully)
pd=read.table("pdata.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

## sanity checks
if(length(unique(pd$folder))>1) {print ("problem: folder name to keep files during pipeline run is not unique")}
outdir=paste(destdir, unique(pd$folder), sep="")
if(length(unique(pd$file1))<dim(pd)[1]) {print ("problem: file1 naming is not unique")}


## if all files are unique, will already rename them now
## else will run prep using original file names and renaming takes place at file combining step (next pipeline)
unames=unique(pd$gsmid)
## check whether multiple files per sample exist
if (length(unames)!=dim(pd)[1]) {fileCombine=T}


## it is ok to have .txt ending in input files (if format is fastq) but for compatibility will exchange this
## this is changed also to the pdata file!
editFileEnd=length(grep(".txt", pd[,1]))>0
if (editFileEnd){
	ufs=unique(pd$file1)
	for (i in 1:length(ufs)){
		system(paste("mv", ufs[i], sub(".txt",".fastq",ufs[i]), sep=" "))		
		}	
	pd[,1]=sub(".txt",".fastq",pd[,1])
	write.table(pd, file="pdata.txt", sep="\t", row.names=F, col.names=T, quote=F)
	}


if (type=="sra"){
	## convert all files listed in pdata from .sra to .fastq format	
	for (i in 1:dim(pd)[1]){
		system(paste("fastq-dump ", cwd, "/", pd$file1[i], " --offset ", offset, sep=""))	
		}
	pd$file1=sub(".sra",".fastq",pd$file1)
		
	}		

## files w/o barcodes
## at this point transfer to outdir where preprocessing continues and create trimconfig files
uf=unique(pd$file1)

if(!fileCombine){
## files are moved to outdir, named based on IDA naming convention (columns 3-10 in pdata file, 11+12 come after alignment) 
## generated trimconfig file just has option to trim

	## make trimconfig file columns
	fm=matrix(nrow=length(unames), ncol=4)
	colnames(fm)=c("file1", "trim_Y_N", "start", "end")

	## go over file by file and add rows to trimconfig and move files to pipeline folder structure
	for (i in 1:length(uf)){
		
		f=pd[i,1]
		ff=paste(pd[i,3],pd[i,4],pd[i,5],pd[i,6],pd[i,7],pd[i,8],pd[i,9],pd[i,10], sep="_")
		if (length(unlist(strsplit(ff, "_")))!=8) {print (paste("problem: IDA naming is not correct for sample", i,sep=" "))}
		fm[i,1]=ff
		fm[i,1]=paste(fm[i,1], ".fastq", sep="")
		if(trim) fm[i,1]=paste(fm[i,1], ".trimmed", sep="")
		if(trima) fm[i,1]=paste(fm[i,1], ".trimmed", sep="")
		fm[i,2]="N"
		fm[i,3]=1
		fm[i,4]=36
		
		
		## move file to pipeline folder structure
		system(paste("mv ", f, " ", outdir, sep=""))		

		## change to file name that matches IDA naming convention
		setwd(outdir)
		system(paste("mv ", f, " ", ff,  ".fastq", sep=""))		
		setwd(cwd)
		}

}


if(fileCombine){
## files are just moved to pipeline folder structure
## trimconfig file generated has option to trim AND option to leave out bad quality runs before combining
## IDA naming compatible final file name is created and written to column 6
	fm=matrix(nrow=dim(pd)[1], ncol=6)
	colnames(fm)=c("file1", "trim_Y_N", "start", "end", "include_Y_N", "IDAsampleName")

	for (i in 1:length(uf)){
		
		f=pd[i,1]
		
		
		fm[i,1]=f
		if (trim) fm[i,1]=paste(fm[i,1], ".trimmed", sep="")
		if(trima) fm[i,1]=paste(fm[i,1], ".trimmed", sep="")
		fm[i,2]="N"
		fm[i,3]=1
		fm[i,4]=36
		fm[i,5]="Y"
		fm[i,6]=paste(pd[i,3],pd[i,4],pd[i,5],pd[i,6],pd[i,7],pd[i,8],pd[i,9],pd[i,10], sep="_")
		if (length(unlist(strsplit(fm[i,6], "_")))!=8) {print (paste("problem: IDA naming is not correct for sample", i,sep=" "))}

		system(paste("mv ", f, " ", outdir, sep=""))		
		

		}

}



## save trimconfig file to destdir/trimconfigs
write.table(fm, file=paste(destdir, "trimconfigs/", unique(pd$folder), "trimconfig.txt" ,sep=""),sep="\t", row.names=FALSE, quote=FALSE)

## save pdata file to outdir
setwd(outdir)
write.table(pd, file="pdata.txt", sep="\t", row.names=F, col.names=T, quote=F)

	
if(editFileEnd){ print("File ending changed to match expected type (fastq)") }
setwd(cwd)
