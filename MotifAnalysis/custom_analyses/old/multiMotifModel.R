

load('/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-06-22/motifScores.acrossTFs.Above0MaxLogodds.withAnno.rda')

#In the meantime, let's get the distal matrix.
distal.f='/srv/gs1/projects/snyder/jzaugg/histoneQTL/hQTL_analysis/results/20140622_newPeaks/DistQTL_list/HiC_local_QTL_dist_QTL.all.txt'
qtl=read.table(distal.f,header=T,sep='\t')
qtl=cbind(qtl,peak1.name=paste(qtl$distQTL,'chr',qtl$chr,':',qtl$peak1.start,'_',qtl$peak1.end,sep=''))

#Annotate peak1 and peak2 - we need this to select the peak from the hmark signal file

first=T
pAnno.t='/srv/gs1/projects/snyder/jzaugg/histoneQTL/peakAnnotation/peakLocations.ALL.HMARK.txt'
signal.t='/srv/gs1/projects/snyder/jzaugg/histoneQTL/hQTL_analysis/PeakData/hMat.norm.ALL.HMARK.peer_lS_5.txt'
hmarks=c('H3K4ME1','H3K4ME3','H3K27AC','RNA','dhs')
for (hmark in hmarks){
    print(hmark)
    curAnno=read.table(gsub('HMARK',hmark,pAnno.t),header=T)
    curAnno=cbind(curAnno[,c('id','chr','start','end')],name=paste(hmark,'chr',curAnno$chr,':',curAnno$start,'_',curAnno$end,sep=''))
    ###### TODO: remove this. very bad.
    curAnno.dupli=which(duplicated(curAnno$name))
    if (length(curAnno.dupli)>0){
       curAnno=curAnno[-curAnno.dupli,]
    }
    if (!first){
       anno=rbind(anno,curAnno)
    }
    if (first){
       anno=curAnno
       first=F
    }   
}
rownames(anno)=anno$name

qtl=cbind(qtl,peak1.id=paste(as.character(qtl$distQTL),'_',as.character(anno[as.character(qtl$peak1.name),'id']),sep=''))
save(qtl,file='/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-06-22/Buffering/HiC_local_QTL_dist_QTL.all.rda')
save(anno,file='/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-06-22/Buffering/combpeakanno.H3K4ME1.H3K4ME3.H3K27AC.dhs.RNA.rda')

signal.t='/srv/gs1/projects/snyder/jzaugg/histoneQTL/hQTL_analysis/PeakData/hMat.norm.ALL.HMARK.peer_lS_5.txt'
hmarks=c('H3K4ME1','H3K4ME3','H3K27AC','RNA')
hmark='H3K4ME1'
curdata=read.table(gsub('HMARK',hmark,signal.t),header=T)
rownames(curdata)=paste(hmark,'_',rownames(curdata),sep='')
peakdata=curdata
hmarks=c('H3K4ME3','H3K27AC','RNA','dhs')
for (hmark in hmarks){
    print(hmark)
    curdata=read.table(gsub('HMARK',hmark,signal.t),header=T)
    rownames(curdata)=paste(hmark,'_',rownames(curdata),sep='')
    all_people=union(colnames(peakdata),colnames(curdata))
    #extend peak data to contain all people
    missed_people=setdiff(all_people,colnames(peakdata))
    peakdata_aug=peakdata
    if (length(missed_people)>0){
       peakdata_aug=data.frame(array(NA,dim=c(nrow(peakdata),length(all_people))))
       rownames(peakdata_aug)=rownames(peakdata)
       colnames(peakdata_aug)=all_people
       peakdata_aug[rownames(peakdata),colnames(peakdata)]=peakdata
    }
    curdata_aug=data.frame(array(NA,dim=c(nrow(curdata),length(all_people))))
    rownames(curdata_aug)=rownames(curdata)
    colnames(curdata_aug)=all_people
    curdata_aug[rownames(curdata),colnames(curdata)]=curdata
    peakdata=rbind(peakdata_aug,curdata_aug)
}
save(peakdata,file='/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-06-22/Buffering/combinedPeakSignal.H3K4ME1.H3K4ME3.H3K27AC.dhs.RNA.rda')

#Add the local QTLs as well. Make a SNPpos - peak file
snppeak=data.frame(chr=qtl$chr,snppos=qtl$dist.snp.pos,peakid=qtl$peak1.id)
snppeak.dupli=which(duplicated(snppeak))
if (length(snppeak.dupli)>0){
   snppeak=snppeak[-snppeak.dupli,]
}
hmarks=c('H3K4ME1','H3K4ME3','H3K27AC','RNA')
for (hmark in hmarks){
    locals=read.table(gsub('HMARK',hmark,'/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-07-28/SNPQTLmatrix.HMARK.maxQTL.gz'),header=TRUE)
    snppeak=rbind(snppeak,data.frame(chr=locals$chr,snppos=locals$snp.position,peakid=paste(hmark,'_',as.character(locals$gene),sep='')))
}
    

#Get all motifs for snp
mScores=cbind(mScores,snp_chrPos=paste(mScores$chr,'_',mScores$SNPpos,sep=''))

##### JOINT MODEL VS SIMPLE MODEL
#==================================
#=====================================
#========================================

#For now, don't care that they're on same SNP
joint_vs_simpleModel=function(y,xdata){
  fitSingleModel=function(idx,y1,xdata1){
    data=rbind(xdata,y1)
    single_model=lm(as.numeric(data[nrow(data),])~as.numeric(data[idx,]))
    return(logLik(single_model))
  }
  
  #Get likelihoods for each single model
  data=rbind(xdata,y)
  sLL=unlist(lapply(1:nrow(xdata),fitSingleModel,y1=data[nrow(data),],xdata1=data[,-nrow(data),]))
  bestm=which.max(sLL)
  #Fit joint model
  data=rbind(xdata,y)
  jointModel=lm(as.numeric(data[nrow(data),])~t(data[-nrow(data),]))
  LRtestStat=2*(logLik(jointModel)-sLL[bestm])
  pval_jointModel=1-pchisq(as.numeric(LRtestStat),nrow(xdata)-1)
  #Return Likelihood ratio test results
  to_return=list()
  to_return[['LRTpvalue']]=pval_jointModel
  to_return[['jointModel']]=jointModel
  to_return[['simpleModel']]=lm(as.numeric(data[nrow(data),])~t(data[bestm,]))
  to_return[['bestSingleMotif']]=rownames(data)[bestm]
  return(to_return)
}


pdf('/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-06-22/Buffering/multiMotifModel.pvals.pdf')
peaks=unique(snppeak$peakid)
hmarks=c('H3K4ME1','H3K4ME3','H3K27AC','RNA')
counter=0
for (hmark_name in hmarks){
    print(counter)
    counter=counter+1
    print(hmark_name)
    curpeaks=peaks[which(grepl(hmark_name,peaks))]
    pvals=c()
    for (peak in curpeaks[1:300]){
    	print(counter)
        counter=counter+1
    	rows_=which(as.character(snppeak$peakid)==peak)
    	snp=paste('chr',snppeak[rows_,'chr'],'_',snppeak[rows_,'snppos'],sep='')
    	hmark=peakdata[peak,]
    	ms=mScores[which(as.character(mScores$snp_chrPos) %in% snp),]
    	if (nrow(ms)>0){
       	   motifdata=ms[,-which(colnames(ms) %in% c('chr','SNPpos','snp_chrPos'))]
       	   jointAnalysis=joint_vs_simpleModel(hmark,motifdata)
       	   #print(jointAnalysis[['LRTpvalue']])
       	   #print(jointAnalysis)
	   pvals=c(pvals,jointAnalysis[['LRTpvalue']])
        }
    }
    hist(pvals,breaks=seq(from=0,to=1,by=0.1),main=hmark_name,xlab='Pvalue - LR test (multi-motif different from single motif model)')
}
dev.off()