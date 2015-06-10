
args=commandArgs(trailingOnly=TRUE)
tfdata.f=args[1]
out=args[2]

#tfdata.f='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/SNPQTLmatrix.RNA.gz.overlapTFBS.gz'
#out='/home/oursu/testpdf'
require(reshape2)
tfdata=read.table(tfdata.f)
tfdata=cbind(tfdata,qtl_and_pass=paste(tfdata[,10],tfdata[,12],sep='_'))
counts=t(table(tfdata[,c(20,18)]))
counts=data.frame(counts)
colnames(counts)=c('TF','QTLtype','count')
count_data=dcast(counts,TF~QTLtype)
count_data=data.frame(count_data,total=rowSums(count_data[,-which(colnames(count_data)=='TF')]))
count_data=count_data[order(-count_data$total),]
write.table(count_data,file=out,quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)

