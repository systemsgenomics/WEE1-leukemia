## R --slave --vanilla --file=$SCRIPTS$(echo qc/filter_trim.R) --args $FILES $PROJ $OFFSET $QP $QMIN $DESTDIR
args <- commandArgs(TRUE)
files=args[1]
proj=args[2]
off=args[3]
qp=args[4]
qm=args[5]
destdir=args[6]

setwd(paste(files, proj, sep=""))



conf=read.table(paste(files,"trimconfigs/", proj, "trimconfig.txt", sep=""), header=TRUE, sep="\t", stringsAsFactors=FALSE)
print(conf)
print("checking if previous intermediates exist for this project, will delete if yes")
if(file.exists(paste(destdir,"intermed_res/",proj, sep=""))) system(paste("rm -r ", destdir,"intermed_res/",proj, sep=""))
system(paste("mkdir ", destdir, "intermed_res/",proj, sep=""))

## check whether files to leave out exist (case: combining several files to one)
if(dim(conf)[2]>5){
	del=conf[,5]=="N"
	conf=conf[!del,]	
	}


for ( i in 1:dim(conf)[1] ){
	print(paste("processing", conf[i,1]))
	if (conf[i,2]=="N"){ 
		## no trimming, filter by quality
				
		system(paste("fastq_quality_filter -i ", conf[i,1]," -o ",destdir, "intermed_res/",proj,"/qfilt_", conf[i,1]," -q ", qm, " -p ", qp," -v -Q ", off, sep=""))
		#system(paste("fastqc ",destdir,"intermed_res/",proj, "/qfilt_", conf[i,1]," --outdir=", destdir, "QC", sep=""))
		}
	if (conf[i,2]=="Y"){ 
		## cut (=trim) read ends, filter by quality
		system(paste("fastx_trimmer -i ", conf[i,1]," -o ", destdir, "intermed_res/",proj,"/trimmed_", conf[i,1], " -f ",conf[i,3]," -l ", conf[i,4]," -v -Q ", off, sep=""))		
		system(paste("fastq_quality_filter -i ", destdir, "intermed_res/",proj,"/trimmed_", conf[i,1]," -o ",destdir, "intermed_res/",proj,"/qfilt_", conf[i,1]," -q ", qm, " -p ", qp," -v -Q ", off, sep=""))
		#system(paste("fastqc ",destdir, "intermed_res/",proj,"/qfilt_", conf[i,1]," --outdir=", destdir, "QC", sep=""))
		}
	print(system("date"))



	

	}

	## combine files if needed
		if(dim(conf)[2]>5){
			setwd(paste(destdir, "intermed_res/",proj, sep=""))
			unames=unique(conf[,6])
			conf[,1]=paste("qfilt_", conf[,1], sep="")
			for (i in 1:length(unames)){
				tocomb=paste(conf[conf[,6]%in%unames[i],1], collapse=" ")
				system(paste("cat ", tocomb, "> qfilt_", unames[i], ".fastq", sep=""))
				system(paste("rm ", tocomb, sep=""))
				}
			}
## changed last line from qfilt_proj_unames to lack proj, added rm step to delete unnecessary files
