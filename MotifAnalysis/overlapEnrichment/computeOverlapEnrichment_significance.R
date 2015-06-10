
args=commandArgs(trailingOnly=TRUE)

tf.f=args[1]
qtlpeak.f=args[2]
out=args[3]

compute_TF_enrichment=function(tf.f,qtlpeak.f,out){
tfdata=read.table(tf.f,header=TRUE,sep='\t')
tfdata[is.na(tfdata)]=0
tfs=setdiff(colnames(tfdata),c('chr','start','end','peak'))
qtlpeak=read.table(qtlpeak.f)
rownames(tfdata)=tfdata$peak
qtlpeaks=as.character(qtlpeak[,1])
#restrict to only peaks with at least 1 TFBS
#peaksWithTFBS=which(rowSums(tfdata[,tfs])>0)
#tfdata=tfdata[peaksWithTFBS,]
#qtlpeaks=as.character(qtlpeak[which(as.character(qtlpeak[,1]) %in% as.character(tfdata$peak)),1])
#tfdata=data.frame(tfdata,QTLpeak=FALSE)
#tfdata[qtlpeaks,'QTLpeak']=TRUE


n=dim(tfdata)[1]
result=data.frame(tf=tfs,
		A=NA,
		Overlap=NA,
		B=length(qtlpeaks),
		Total=n,
		enrichment=NA,
		enrichment_p=NA,
		confLow=NA,confHigh=NA)
rownames(result)=tfs
for (tf in tfs){
    tf_peaks=as.character(tfdata[which(as.numeric(as.character(tfdata[,tf]))>0),'peak'])
    qtl_peaks_with_tf=intersect(tf_peaks,qtlpeaks)
    m=data.frame(TF=c(length(qtl_peaks_with_tf),length(tf_peaks)-length(qtl_peaks_with_tf)),
		notTF=c(length(qtlpeaks)-length(qtl_peaks_with_tf),n-length(tf_peaks)-length(qtlpeaks)+length(qtl_peaks_with_tf)))
    ftest=fisher.test(m)
    result[tf,'enrichment']=ftest$estimate
    result[tf,'enrichment_p']=ftest$p.value
    result[tf,'A']=length(tf_peaks)
    result[tf,'Overlap']=length(qtl_peaks_with_tf)
    result[tf,'confLow']=ftest$conf.int[1]
    result[tf,'confHigh']=ftest$conf.int[2]
}
result=result[order(result$enrichment_p),]
result=data.frame(result,BH=p.adjust(result$enrichment_p))
rownames(result)=NULL
write.table(result,
	file=out,
	sep='\t',quote=F,row.names=F,col.names=T)
return(result)
}

### COMPUTE ENRICHMENTS
#######################
print('Computing enrichment')
result=compute_TF_enrichment(tf.f,qtlpeak.f,out)
#replace complicated TF names by TF HGNCs
#hgnc.f='/home/oursu/ENCODE.hg19.TFBS.QC.metadata.jun2012 - TFs_SPP_pooled.tsv'
#hgnc=read.csv(hgnc.f,sep='\t')
#rownames(hgnc)=hgnc$FILENAME
result2=data.frame(result,TFname=result$tf)

### PLOT ENRICHMENT
###################
print('Plotting')
require(ggplot2)
pdf(paste(out,'.pdf',sep=''),height=12,width=10)
res.cur=result2
res.cur$TFname=factor(res.cur$TFname,levels=res.cur$TFname[order(res.cur$enrichment)])
res.cur=cbind(res.cur,sig_after_BHcorrection=(as.numeric(as.character(res.cur$BH))<=0.05))
print(ggplot(res.cur, aes(y=enrichment, x=TFname,col=sig_after_BHcorrection)) + coord_flip()+geom_point()+geom_errorbar(aes(ymax = res.cur$confLow, 
ymin=res.cur$confHigh))+ylim(0,4)+
ggtitle(basename(out)))
print(ggplot(res.cur, aes(y=enrichment, x=TFname,col=sig_after_BHcorrection)) + coord_flip()+geom_point()+geom_errorbar(aes(ymax = res.cur$confLow, 
ymin=res.cur$confHigh))+ylim(0,20)+
ggtitle(basename(out)))
dev.off()



