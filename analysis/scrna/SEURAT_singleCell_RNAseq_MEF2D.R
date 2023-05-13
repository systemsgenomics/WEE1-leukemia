## --- Wrapper script to run Seurat analysis



## conda activate R4_hdf5_seurat4
library(Seurat)
library(Matrix)


##

readRaw=function(inputDir){
    
    matrix_dir = inputDir
    barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
    features.path <- paste0(matrix_dir, "features.tsv.gz")
    matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
    mat <- readMM(file = matrix.path)
    feature.names = read.delim(features.path,
    header = FALSE,
    stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path,
    header = FALSE,
    stringsAsFactors = FALSE)
    colnames(mat) = barcode.names$V1
    ## matching to symbols not unique, fix here
    feature.names[duplicated(feature.names[,2]),2]=paste(feature.names[duplicated(feature.names[,2]),2],feature.names[duplicated(feature.names[,2]),1], sep=":")
    
    #rownames(mat) = feature.names$V1
    rownames(mat) = feature.names$V2
    
    dat.r=mat[feature.names[,3]=="Gene Expression",]
    s <- CreateSeuratObject(dat.r)

    if(nrow(dat.r)<length(feature.names[,3])){
        isHashAb=1:nrow(mat)%in%grep("Hash",rownames(mat))
        dat.p=mat[!isHashAb&feature.names[,3]=="Antibody Capture",]
        s[["ADT"]] <- CreateAssayObject(counts = as.matrix(dat.p))

        if(any(isHashAb)){
            dat.h=mat[isHashAb,]
            s[["HTO"]] <- CreateAssayObject(counts = as.matrix(dat.h))
        }
    }
    return(s)
    
}


## --- PREPARE INPUT data


# start from raw cellranger output
## diag BM
diag="cellranger/filt_mtx/L5/10k/"
## day2 blood
d2="cellranger/filt_mtx/L10/filtered_feature_bc_matrix/"


dd1=readRaw(diag)
dd2=readRaw(d2)

## metadata, diag
m1=read.table("adata1_clu_0.csv", sep=",", header=T)
rownames(m1)=m1[,1]
colnames(m1)[2]="scanpy_leidenr1_cluster"
## leukemic cell clusters
m1=m1[m1[,2]%in%c(0,1,2,3,4,5,6),]

## metadata, day2
m2=read.table("adata2_clu_0.csv", sep=",", header=T)
rownames(m2)=m2[,1]
colnames(m2)[2]="scanpy_leidenr1_cluster"
## leukemic cell clusters and those expressing CD19 or CD20
m2.l=m2[m2[,2]%in%c(0,1,2,3,4,5,6),]


## subset to leukemic cells and other B-lineage cells in data
dd1=dd1[,colnames(dd1)%in%m1[,1]]
dd2.l=dd2[,colnames(dd2)%in%m2.l[,1]]


dd1=AddMetaData(dd1, m1)
dd2.l=AddMetaData(dd2.l, m2.l)

dd1=AddMetaData(dd1, rep("diag", ncol(dd1)),col.name="timepoint")
dd2.l=AddMetaData(dd2.l, rep("day2", ncol(dd2.l)),col.name="timepoint")


## Create list of separate Seurat objects ready for normalization ----

N.list=list(dd1,dd2.l)
names(N.list)=c("10X_diag","10X_d2_leukemic")

## LogNormalize
N.list.LogNorm <- lapply(X = N.list, FUN = function(x) {
    #    x <- x[match(keep_genes, rownames(x)), ]
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
    x <- ScaleData(x)
    x <-RunPCA(x, ndims.print=1, nfeatures.print=1)
    
})

pid=names(N.list)
for(i in 1:length(pid)) saveRDS(N.list.LogNorm[[i]], file=paste(pid[i],"LogNorm.rds", sep="_"))
