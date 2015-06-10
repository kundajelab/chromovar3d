args=commandArgs(trailingOnly=TRUE)


motifscore.f='/srv/gsfs0/projects/kundaje/users/oursu/histoneQTLproject/motif_analysis/data/motifPWM/matchesPerTF/correctedByPouyaControls/motifFreqMatrices/matchesPerTF/motif.pouya.Motif.MOTIF_known1scanThresh2.MATCH.logOdds2.gz'
controlscore.f='/srv/gsfs0/projects/kundaje/users/oursu/histoneQTLproject/motif_analysis/data/motifPWM/matchesPerTF/correctedByPouyaControls/motifFreqMatrices/matchesPerTF/motif.pouya.Motif.MOTIF_known1scanThresh2.MATCH.logOdds2.gzControlScores.gz'

out='/srv/gsfs0/projects/kundaje/users/oursu/histoneQTLproject/motif_analysis/data/motifPWM/matchesPerTF/correctedByPouyaControls/plots/_'

tf_names='SPI1,BRCA1,NR4A,TCF12,TCF7L2,GATA,HEY1,HSF,MXI1,NFKB,NRF1,SP1,TATA,ZBTB33,STAT,ZBTB7A,AP1,NFE2,NFY,SREBP,EGR1,ELF1,NANOG,PBX3,TFAP2,IRF,ZEB1,ATF3,BATF,BHLHE40,CEBPB,EP300,MYC,NR2C2,PRDM1,E2F,FOXA,POU5F1,ETS,MAF,MEF2,TAL1,RFX5,SIX5,YY1,SRF,CTCF,ESRRA,HNF4,NR3C1,POU2F2,ZNF143,REST,RXRA,EBF1,PAX5'
tfs=strsplit(tf_names,',')[[1]]

#MOTIFSCORES
pdf(paste(out,'motifScoreDensity.pdf',sep=''))
motifScores_v=list()
motifScores=list()
for (tf in tfs){
    print(tf)
    data=read.table(gsub('MOTIF',tf,motifscore.f))
    motifScores[[tf]]=ecdf(data[,5])	
    motifScores_v[[tf]]=data[,5]
    plot(density(data[,5]),xlim=c(2,20),ylim=c(0,1),col='gray',main='',xlab='Homer Motif score')
    par(new=T)
}
dev.off()

pdf(paste(out,'motifScoreECDF.pdf',sep=''))
colorcounter=1
colorlabels=c()
for (tf in names(motifScores)){
    print(tf)
    plot(motifScores[[tf]],xlim=c(2,20),ylim=c(0,1),col=colorcounter,main='',xlab='Homer Motif score')
    par(new=T)
    colorcounter=colorcounter+1
    colorlabels=c(colorlabels,tf)
}
colorcounter=colorcounter-1
#legend("topleft",fill=c(1:colorcounter),legend=colorlabels)
dev.off()
save(motifScores,file=paste(out,'MotifScores.rda',sep=''))
save(motifScores_v,file=paste(out,'MotifScores_value.rda',sep=''))

#CORRECTED MOTIFSCORES
controlScores=list()
controlScores_v=list()
maxControl=list()
#pdf(paste(out,'motifScoreMinusControlDensity.pdf',sep=''))
for (tf in tfs){
    print(tf)
    data=read.table(gsub('MOTIF',tf,controlscore.f))
    data=cbind(data,maxVal=apply(data[,c(6:ncol(data))],1,max),featureName=gsub('-2$','',data[,1]))
    data=data[order(data$maxVal,decreasing=T),]
    dupli_data=which(duplicated(data$featureName))
    if (length(dupli_data)>0){
       data=data[-which(duplicated(data$featureName)),]
    }
    rownames(data)=data$featureName
    true_data=read.table(gsub('MOTIF',tf,motifscore.f))
    #true data contains dubly reported motifs, so we can just take 1 version
    dupli_true=which(duplicated(paste(true_data[,1],'_',true_data[,2],'_',true_data[,3],sep='')))	
    if (length(dupli_true)>0){
       true_data=true_data[order(true_data[,5],decreasing=T),]
       true_data=true_data[-which(duplicated(paste(true_data[,1],'_',true_data[,2],'_',true_data[,3],sep=''))),]
    }
    rownames(true_data)=paste(true_data[,1],':',true_data[,2],'-',true_data[,3],sep='')
    common_rows=intersect(rownames(data),rownames(true_data))
    controlScores[[tf]]=ecdf(true_data[common_rows,5]-data[common_rows,'maxVal'])
    #controlScores_v[[tf]]=true_data[common_rows,5]-data[common_rows,'maxVal']
    maxControl[[tf]]=ecdf(data[common_rows,'maxVal'])
    #plot(density(controlScores_v[[tf]]),xlim=c(2,20),ylim=c(0,1),col='gray',main='',xlab='Homer Motif score - max control score')
    #par(new=T)
}
#dev.off()
pdf(paste(out,'controlScoreECDF.pdf',sep=''))
colorcounter=1
colorlabels=c()
for (tf in names(controlScores)){
    print(tf)
    plot(controlScores[[tf]],xlim=c(-10,30),ylim=c(0,1),col=colorcounter,main='',xlab='Homer Motif score - max control score at the site')
    par(new=T)
    colorcounter=colorcounter+1
    colorlabels=c(colorlabels,tf)
}
dev.off()

pdf(paste(out,'maxcontrolScoreECDF.pdf',sep=''))
colorcounter=1
colorlabels=c()
for (tf in names(controlScores)){
    print(tf)
    plot(maxControl[[tf]],xlim=c(-20,30),ylim=c(0,1),col=colorcounter,main='',xlab='Max control score at the site')
    par(new=T)
    colorcounter=colorcounter+1
    colorlabels=c(colorlabels,tf)
}
dev.off()


pdf(paste(out,'maxcontrolScoreECDF.pdf',sep=''))
colorcounter=1
colorlabels=c()
for (tf in names(controlScores)){
    print(tf)
    plot(maxControl[[tf]],xlim=c(-20,30),ylim=c(0,1),col=colorcounter,main='',xlab='Max control score at the site')
    par(new=T)
    colorcounter=colorcounter+1
    colorlabels=c(colorlabels,tf)
}
dev.off()

colorcounter=colorcounter-1
#legend("topleft",fill=c(1:colorcounter),legend=colorlabels)
dev.off()
save(controlScores,file=paste(out,'MotifMinusMaxControlScores.rda',sep=''))
save(maxControl,file=paste(out,'MaxControlScores.rda',sep=''))

#Make some files with the motif thresholds for each TF mased on motif score alone
ms=c(0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5)
tf_score=data.frame(TF=tfs,scoreThreshQ0=2)
rownames(tf_score)=tf_score$TF
for (m in ms){
    cur_values=c()
    for (i in c(1:nrow(tf_score))){
    	cur_values=c(cur_values,quantile(motifScores[[tfs[i]]],m))
    }
    tf_score=cbind(tf_score,cur_values)
    colnames(tf_score)[ncol(tf_score)]=paste('scoreThreshQ',as.character(m),sep='')
}
write.table(tf_score,file='/srv/gsfs0/projects/kundaje/users/oursu/histoneQTLproject/motif_analysis/data/motifPWM/matchesPerTF/correctedByPouyaControls/MatchThresholds/MotifQuantileThresholds.SITES_motifScanLogOdds2onRefGenome.notControlCorrected.tab',sep='\t',quote=F,row.names=F,col.names=T)
require(pheatmap)
pdf('/srv/gsfs0/projects/kundaje/users/oursu/histoneQTLproject/motif_analysis/data/motifPWM/matchesPerTF/correctedByPouyaControls/MatchThresholds/MotifQuantileThresholds.SITES_motifScanLogOdds2onRefGenome.notControlCorrected.pdf')
pheatmap(as.matrix(tf_score[,-1]),cluster_cols=F,cellheight=1,cellwidth=1,fontsize=1)
dev.off()













