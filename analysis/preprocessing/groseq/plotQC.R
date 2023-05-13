d=read.table("tagCountDistribution.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
d=as.matrix(d)
d[,2]=d[,2]*100

if(dim(d)[1]>2){
	pdf(file="tagCountDistribution.pdf")
	barplot(d[2:dim(d)[1],2], main="tag count distribution", xlab="position", ylab="% tags", names.arg=1:(dim(d)[1]-1))
	dev.off()
	}
## else no plotting because tag count has been set to 1

d=read.table("tagAutocorrelation.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
d=as.matrix(d)
#d=d[999:3001,]
d=d[1501:2501,]
colnames(d)=c("distance in bp", "same strand", "opposite strand")
same=cbind(d[,1],d[,2])
opp=cbind(d[,1],d[,3])
plot_colors <- c(rgb(r=0.0,g=0.0,b=0.9), "red", "forestgreen",rgb(r=.5,g=0,b=0.5))
pdf(file="tagAutocorrelation.pdf")
plot(same, type="l", col=plot_colors[1], ylim=range(d[,2:3]), xlab="Relative distance between reads (bp)", ylab="Total read pairs", cex.lab=0.8, lwd=2)
lines(opp, type="l", lty=2, lwd=2, col=plot_colors[2])
legend("bottomright", c("same strand", "opposite strand"), cex=0.8, col=plot_colors, lty=1:2, lwd=2, bty="n")
dev.off()

d=read.table("tagGCcontent.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
d=as.matrix(d)
d2p=cbind(d[,1],d[,3])
g=read.table("genomeGCcontent.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
g=as.matrix(g)
g2p=cbind(g[,1],g[,3])
pdf(file="tagGCcontent.pdf")
plot(d2p, type="l", col=plot_colors[1],  xlab="GC-content of fragments", ylab="Normalized fraction", cex.lab=0.8, lwd=2)
lines(g2p, type="l", lty=2, lwd=2, col=plot_colors[2])
legend("topright", c("sample","genome"), cex=0.8, col=plot_colors, lty=1:2, lwd=2, bty="n")
dev.off()

d=read.table("tagFreqUniq.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
d=as.matrix(d)
A=cbind(d[,1],d[,2])
C=cbind(d[,1],d[,3])
G=cbind(d[,1],d[,4])
T=cbind(d[,1],d[,5])

pdf(file="tagFreqUnique.pdf")
plot(A, type="l", col=plot_colors[1], ylim=range(d[,2:3]),  xlab="Distance from 5'end of reads (bp)", ylab="Nucleotide frequency", cex.lab=0.8, lwd=2)
lines(C, type="l", lty=2, lwd=2, col=plot_colors[2])
lines(G, type="l", lty=2, lwd=2, col=plot_colors[3])
lines(T, type="l", lty=2, lwd=2, col=plot_colors[4])
legend("bottomright", colnames(d[,2:5]), cex=0.8, col=plot_colors, lty=1:4, lwd=2, bty="n")
dev.off()

d=read.table("tagFreq.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
d=as.matrix(d)
A=cbind(d[,1],d[,2])
C=cbind(d[,1],d[,3])
G=cbind(d[,1],d[,4])
T=cbind(d[,1],d[,5])

pdf(file="tagFreq.pdf")
plot(A, type="l", col=plot_colors[1], ylim=range(d[,2:3]),  xlab="Distance from 5'end of reads (bp)", ylab="Nucleotide frequency", cex.lab=0.8, lwd=2)
lines(C, type="l", lty=2, lwd=2, col=plot_colors[2])
lines(G, type="l", lty=2, lwd=2, col=plot_colors[3])
lines(T, type="l", lty=2, lwd=2, col=plot_colors[4])
legend("bottomright", colnames(d[,2:5]), cex=0.8, col=plot_colors, lty=1:4, lwd=2, bty="n")
dev.off()
