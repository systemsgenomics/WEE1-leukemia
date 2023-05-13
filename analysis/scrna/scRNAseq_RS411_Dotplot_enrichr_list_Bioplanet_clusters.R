# Dotplot visualization
#scRNAseq RS411 enrichr from Bioplanet sorted by pval

library(ggplot2)

## can drop background grid with this 
library(cowplot)
theme_cowplot()



#take top5 each, uniq
eterms=unique(system("head -5 top5/BioPlanet_2019*mod | cut -f 1 | uniq ", intern=T))
####get weird names and empty rows out
del=grep("==>",eterms)
eterms=eterms[-del]
eterms=eterms[!eterms==""]

#eterms=eterms[-1]
ENC=system("ls BioPlanet_2019*cleaned", intern=T)
etable=as.data.frame(matrix(nrow=length(eterms), ncol=3*length(ENC)), stringsAsFactors=F)

rownames(etable)=eterms
for (i in 1:length(ENC)){
  e=read.table(ENC[i], header=T, sep="\t", stringsAsFactors=F)
  e=e[,c(1,3,4,2)]
  e[,4]=sub("\\/.*","",e[,4])
  e=e[e$Term%in%eterms,]
  ematch=match(rownames(etable),e[,1])
  etable[,(3*(i-1)+1):(3*i)]=e[ematch,2:4]
  colnames(etable)[(3*(i-1)+1):(3*i)]=paste(sub(".*BioPlanet_2019_table_","",ENC[i]),colnames(e)[2:4],sep=";")
}
#this table has all the results. supplement table. dotplot is from cutoff applied

write.table(cbind(eterms,etable),file="temp5BioPlanet_overlap_table_scRNAseq_RS411_090620.txt", row.names=F, col.names=T, quote=F, sep="\t")
ap=grep("Adjusted",colnames(etable))

res=etable
rownames(res)=paste(1:nrow(res),rownames(res),sep="_")
features.plot=rownames(res)
data.plot=res[,1:3]
colnames(data.plot)=sub(".*;","",colnames(data.plot))
id=rep(sub(".txt.*","",colnames(res)[1]),nrow(res))
for(i in 2:length(ENC)){
  tt=res[,(3*(i-1)+1):(3*i)]
  id=c(id,rep(sub(".txt.*","",colnames(res)[(3*i)]),nrow(res)))
  colnames(tt)=colnames(data.plot)
  data.plot=rbind(data.plot,tt)
}

data.plot=cbind(rep(features.plot,length(ENC)),data.plot,id)
colnames(data.plot)[1]="Terms"

#cut scale #add scale here if pvalue goes smaller than --> same red color
#pval[pval<1e-9]=1e-9
#more significant darker color in picture
pval=-10*log10(data.plot$Adjusted.P.value)
pval[is.na(pval)]=1 #jotku on 0 p-val
summary(pval)

pval[pval>20]=20



#make continous value from ovelapped genes
data.plot$Overlap=as.numeric(data.plot$Overlap)

#if not any gene as overlap. adjust to small value for viz purposes
data.plot$Overlap[is.na(data.plot$Overlap)]=0.0001

#alphabetic order
data.plot$Terms=factor(as.vector(data.plot$Terms), levels=rev(unique(data.plot$Terms)))

pdf("BioPlanet_pathways_scRNAseq_RS411_uniq_gene_names_top5_090620.pdf", height = 10)
plot <- ggplot(data = data.plot, mapping = aes_string(x = "id", y = "Terms"))
plot <- plot + geom_point(mapping = aes_string(size = "Overlap",colour=pval))
plot <- plot + scale_size(range = c(1, 10))
#plot <- plot + geom_point(aes(colour=pval))
plot <- plot + scale_color_gradient(low="grey", high="red")
print(plot)
dev.off()


