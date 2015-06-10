
fileregex='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/pairAnalysis/TF_TF/*summary.txt'
outf='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/pairAnalysis/TF_TF/ZscoreHeatmap.pdf'
fs=Sys.glob(fileregex)
THRESH=10

first=TRUE
for (f in fs){
	curdata=read.table(f,header=TRUE)
	#curdata[,'Total']=curdata$Total+1
	#curdata[,'F1atsite1']=curdata$F1atsite1+1
	#curdata[,'F2atsite2']=curdata$F2atsite2+1
	#curdata[,'Paired1atsite1with2atsite2']=curdata$Paired1atsite1with2atsite2+1
	if (first!=TRUE){
		data=rbind(data,curdata)
	}
	if (first==TRUE){
		data=curdata
		first=FALSE
	}
}
data[intersect(which(data$Item1=='MergedPeaks_ChIPseq_ZZZ3'),which(data$Item2=='MergedPeaks_ChIPseq_ZNF274')),]
#compute parameters for the zscore, and the zscore
data=data.frame(data,
	m12=data$F1atsite1*data$F2atsite2/(data$Total^2),
	m21=data$F1atsite2*data$F2atsite1/(data$Total^2),
	sd12=((1/data$Total)*(data$F1atsite1*data$F2atsite2/(data$Total^2))*(1-data$F1atsite1*data$F2atsite2/(data$Total^2)))^0.5,
	sd21=((1/data$Total)*(data$F1atsite2*data$F2atsite1/(data$Total^2))*(1-data$F1atsite2*data$F2atsite1/(data$Total^2)))^0.5)
data=data.frame(data,
	z12=(data$Paired1atsite1with2atsite2/data$Total-data$m12)/data$sd12,
	z21=(data$Paired1atsite2with2atsite1/data$Total-data$m21)/data$sd21)
#filter out things with <10 at either side
filter_out=union(which(data$Paired1atsite1with2atsite2<=THRESH),union(which(data$F2atsite2<=THRESH),
	which(data$F1atsite1<=THRESH)))
data[filter_out,'z12']=0
require(reshape2)
itemxitem=dcast(data, Item1 ~ Item2, value.var="z12")
rownames(itemxitem)=itemxitem[,1]
itemxitem=itemxitem[,-1]
itemxitem[abs(itemxitem)<=1.65]=0

#heatmap the zscores
pdf(outf)
require(pheatmap)
pheatmap(itemxitem[1:70,1:70],cellheight=1,cellwidth=1,fontsize=1)
dev.off()

#heatmap the sig zscores
