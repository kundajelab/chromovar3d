
args=commandArgs(trailingOnly=TRUE)
require(GenomicRanges)
require(rtracklayer)


motifScores.f=args[1]
tfbs.tf.f=args[2] #chip-seq
tfbs.total.f=args[3] #chip-seq all

#motifScores.f='/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/data/motifPWM/matchesPerTF/correctedByPouyaControls/motifFreqMatrices/matchesPerTF/motif.pouya.Motif.SP1_known1scanThresh2.MATCH.logOdds2.gz'
#tfbs.tf.f='/srv/gsfs0/projects/kundaje/commonRepository/encode/data/byDataType/peaks_spp/mar2012/distinct/idrOptimalBlackListFilt/wgEncodeHaibTfbsGm12878Sp1Pcr1xAlnRep0.bam_VS_wgEncodeHaibTfbsGm12878RxlchPcr1xAlnRep0.bam.regionPeak.gz'
#tfbs.total.f='/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/data/TFBS_GM12878.bed.gz'

motifScores=read.table(motifScores.f)
motifScores.dupli=which(duplicated(paste(motifScores[,1],'_',motifScores[,2],'_',motifScores[,3],sep='')))
if (length(motifScores.dupli)>0){
   motifScores=motifScores[-motifScores.dupli,]
}
controlScores=read.table(paste(motifScores.f,'ControlScores.gz',sep=''))

### quick plot of observed vs control scores
pdf('/srv/gsfs0/projects/kundaje/users/oursu/testAlignments/testMotifFDR_CTCF.pdf')
par(mfrow=c(2,1))
hist(motifScores[,5],col='blue',xlim=c(-10,20),breaks=8)
hist(as.vector(as.matrix(controlScores[,c(6:ncol(controlScores))])),col='pink',xlim=c(-10,20),breaks=40)
dev.off()

allm=as.vector(as.matrix(controlScores[,c(6:ncol(controlScores))]))
m=motifScores[,5]
for (i in c(-10:20)){
    print(i)
    controls=ncol(controlScores)-5
    print((length(which(allm>i))/controls)/(length(which(m>i))))
}

controlScores.dupli=which(duplicated(controlScores[,1]))
if (length(controlScores.dupli)>0){
   controlScores=controlScores[-controlScores.dupli,]
}
rownames(controlScores)=controlScores[,1]
controlScores=controlScores[paste(motifScores[,1],':',motifScores[,2],'-',motifScores[,3],sep=''),]
controlScore_cols=c(6:ncol(controlScores))
maxControlScores=apply(controlScores[,controlScore_cols],1,max)

motifCol=5
tfbs.tf=read.table(tfbs.tf.f)
tfbs.total=read.table(tfbs.total.f)


motifScores.gr=GRanges(seqnames =Rle(as.character(motifScores[,1])), ranges =IRanges(motifScores[,2], end = motifScores[,3], 
				names = paste(motifScores[,1],'_',motifScores[,2],'_',motifScores[,3],sep='')),score = motifScores[,5])
controlScores.gr=GRanges(seqnames =Rle(as.character(motifScores[,1])), ranges =IRanges(motifScores[,2], end = motifScores[,3],
                                names = paste(motifScores[,1],'_',motifScores[,2],'_',motifScores[,3],sep='')),score = maxControlScores)
tfbs.tf.gr=GRanges(seqnames =Rle(as.character(tfbs.tf[,1])), ranges =IRanges(tfbs.tf[,2], end = tfbs.tf[,3],
                                names = paste(tfbs.tf[,1],'_',tfbs.tf[,2],'_',tfbs.tf[,3],sep='')))
tfbs.total.gr=GRanges(seqnames =Rle(as.character(tfbs.total[,1])), ranges =IRanges(tfbs.total[,2], end = tfbs.total[,3],
                                names = paste(tfbs.total[,1],'_',tfbs.total[,2],'_',tfbs.total[,3],sep='')))

# ------  TEST: Check thar tfbs.tf are completely contained inside tfbs.total
tfbs.total.OVERLAP.tf=as.data.frame(findOverlaps(query=tfbs.tf.gr,subject=tfbs.total.gr,select='all',type='within'))
colnames(tfbs.total.OVERLAP.tf)=c('TF','allTFBS')
if (dim(tfbs.total.OVERLAP.tf)[1]<dim(tfbs.tf)[1]){
   print("TF TFBS are not fully contained in the merged TFBS")
}
#-------

#Label TFBS by whether they have the TF or not.
tfbs.total=cbind(tfbs.total,TFpeak=FALSE)
tfbs.total[tfbs.total.OVERLAP.tf[,'allTFBS'],'TFpeak']=TRUE

#For each TFBS, get the maximum motif score
tfbs.total.OVERLAP.motifs=as.data.frame(findOverlaps(query=motifScores.gr,subject=tfbs.total.gr,select='all',type='within'))
tfbs.total.OVERLAP.motifs=cbind(tfbs.total.OVERLAP.motifs,motifScore=as.data.frame(motifScores.gr)[,'score'],maxControlScore=as.data.frame(controlScores.gr)[,'score'])
tfbs.total.OVERLAP.motifs=tfbs.total.OVERLAP.motifs[order(tfbs.total.OVERLAP.motifs,decreasing=TRUE),]
tfbs.total.OVERLAP.motifs.dupli=which(duplicated(tfbs.total.OVERLAP.motifs[,2]))
tfbs.total.OVERLAP.motifs.dupli=c(tfbs.total.OVERLAP.motifs.dupli,which(is.na(tfbs.total.OVERLAP.motifs[,2])))
if (length(tfbs.total.OVERLAP.motifs.dupli)>0){
   tfbs.total.OVERLAP.motifs=tfbs.total.OVERLAP.motifs[-tfbs.total.OVERLAP.motifs.dupli,]
}

tfbs.total=cbind(tfbs.total,maxMotifScore=NA,maxControlScore=NA)
rownames(tfbs.total.OVERLAP.motifs)=tfbs.total.OVERLAP.motifs[,2]
common=intersect(rownames(tfbs.total.OVERLAP.motifs),rownames(tfbs.total))
tfbs.total[tfbs.total.OVERLAP.motifs[common,2],'maxMotifScore']=motifScores[tfbs.total.OVERLAP.motifs[common,1],motifCol]
tfbs.total[tfbs.total.OVERLAP.motifs[common,2],'maxControlScore']=maxControlScores[as.numeric(common)]

#Now, we remove peaks without scores
notScored=which(is.na(tfbs.total$maxMotifScore))
tfbs.total.scored=tfbs.total
if (length(notScored)>0){
   tfbs.total.scored=tfbs.total[-notScored,]
}

#Get all peaks and same number of non-peaks
numExamples=table(tfbs.total.scored$TFpeak)
minNumber=min(numExamples)
print(paste("making example sets of size ",minNumber,sep=''))
tfbs.total.scored.shuffled=tfbs.total.scored[sample(nrow(tfbs.total.scored)),]
examples=tfbs.total.scored.shuffled[union(which(tfbs.total.scored.shuffled$TFpeak==TRUE)[1:minNumber],which(tfbs.total.scored.shuffled$TFpeak==FALSE)[1:minNumber]),]
examples=cbind(examples,class=1*examples$TFpeak)

#Finally, logistic regression!
TFglm=glm(class ~ maxMotifScore, data = examples, family = binomial(logit))
pdf(paste(motifScores.f,'LogisticReg.plot.pdf',sep=''))
plot(examples$maxMotifScore,examples$class,xlab="Max motif score in the peak",ylab="Probability of being a peak for TF") 
curve(predict(TFglm,data.frame(maxMotifScore=x),type="resp"),add=TRUE) # draws a curve based on prediction from logistic regression model
plot(density(examples[which(examples$class==0),'maxMotifScore']),xlim=c(min(examples$maxMotifScore),max(examples$maxMotifScore)),ylim=c(0,3),xlab='',ylab='',main='')
par(new=T)
plot(density(examples[which(examples$class==1),'maxMotifScore']),xlim=c(min(examples$maxMotifScore),max(examples$maxMotifScore)),ylim=c(0,3),col='red',
xlab="Max motif score in the ChIP-seq peak",ylab="density",main="Distribution of motif scores in ChIP-seq peaks\n for TF peaks vs negative set (other peaks)")
legend('topleft',fill=c('black','red'),legend=c('Negatives','Positives'))
hist(examples[which(examples$class==1),'maxMotifScore'], col=rgb(1,0,0,0.5),xlim=c(0,max(examples$maxMotifScore)), main="Distribution of motif scores in ChIP-seq peaks\n for TF peaks vs negative set (other peaks)",xlab="Max motif score in the peak")
hist(examples[which(examples$class==0),'maxMotifScore'], col=rgb(0,0,1,0.5), add=T) #breaks=seq(from=2,to=max(examples$maxMotifScore),by=0.5),add=T)
legend('topleft',fill=c('blue','red'),legend=c('Negatives','Positives'))
box()
require(popbio)
logi.hist.plot(examples$maxMotifScore,examples$class,boxp=FALSE,type="hist",col="gray")
dev.off()

TFglm_withControl=glm(class ~ maxMotifScore+maxControlScore, data = examples, family = binomial(logit))
pdf(paste(motifScores.f,'LogisticReg.plot.withControl.pdf',sep=''))
xlabel=paste(as.character(TFglm_withControl$coefficients[['maxMotifScore']]),'*maxMotifScore + ',as.character(TFglm_withControl$coefficients[['maxControlScore']]),'*maxControlScore',sep='')
values_=TFglm_withControl$coefficients[['maxMotifScore']]*examples$maxMotifScore+TFglm_withControl$coefficients[['maxControlScore']]*examples$maxControlScore
plot(values_,examples$class,xlab=xlabel,ylab='p(bound)')
require(popbio)
logi.hist.plot(values_,examples$class,boxp=FALSE,type="hist",col="gray")
hist(values_[which(examples$class==1)], col=rgb(1,0,0,0.5),xlim=c(0,max(values_)), main="Distribution of motif scores in ChIP-seq peaks\n for TF peaks vs negative set (other peaks)",xlab=xlabel)
hist(values_[which(examples$class==0)], col=rgb(0,0,1,0.5), add=T) #breaks=seq(from=2,to=max(examples$maxMotifScore),by=0.5),add=T)
dev.off()


#Finally, find out the motif threshold.
motifThresh=data.frame(motifThresh=NA,motifThresh_withControl=NA)
motifThresh[1,'motifThresh']=(-TFglm$coefficients[['(Intercept)']])/TFglm$coefficients[['maxMotifScore']]
motifThresh[1,'motifThresh_withControl']=-TFglm_withControl$coefficients[['(Intercept)']]
write.table(motifThresh,file=paste(motifScores.f,'LogisticReg.threshold',sep=''),quote=F,sep='\t',row.names=F,col.names=T)

