
mScores.t='/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-06-22/motifScores/pouya.PWM.motif.Homer.Motif.theTF_known1/pouya.PWM.motif.Homer.Motif.theTF_known1INsnp_scores.HMARK_motifLen*.ALL.TFBS_GM12878.bed.gz.motif.score.avg.rda'
tfs.text='SPI1,BRCA1,NR4A,TCF12,TCF7L2,GATA,HEY1,HSF,MXI1,NFKB,NRF1,SP1,TATA,ZBTB33,STAT,ZBTB7A,AP1,NFE2,NFY,SREBP,EGR1,ELF1,NANOG,PBX3,TFAP2,IRF,ZEB1,ATF3,BATF,BHLHE40,CEBPB,EP300,MYC,NR2C2,PRDM1,E2F,FOXA,POU5F1,ETS,MAF,MEF2,TAL1,RFX5,SIX5,YY1,SRF,CTCF,ESRRA,HNF4,NR3C1,POU2F2,ZNF143,REST,RXRA,EBF1,PAX5'
tfs=strsplit(tfs.text,',')[[1]]
tfs=setdiff(tfs,'EP300')

#Try to put all motif scores into 1 giant matrix. 
first=T
for (tf in tfs){
    print(tf)
    load(Sys.glob(gsub('theTF',tf,gsub('HMARK','H3K4ME3',mScores.t))))	
    cur_mScore=motif.score.avg
    motifMax=apply(cur_mScore,1,max)
    cur_mScore=cur_mScore[which(motifMax>=0),]
    rownames(cur_mScore)=paste(tf,'.',rownames(cur_mScore),sep='')
    if (first){
       mScores=cur_mScore
       first=F
    }
    if (first==F){
       mScores=rbind(mScores,cur_mScore)
    }   
}
#If this works, we'll save this motif matrix for easy access later on.
save(mScores,file='/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-06-22/motifScores.acrossTFs.Above0MaxLogodds.rda')

#Add in SNP chr, position and TF

chromos=unlist(regmatches(rownames(mScores),gregexpr("chr[0-9]*",rownames(mScores))))
mScores2=cbind(mScores,chr=chromos)
starts=gsub('_','',unlist(regmatches(rownames(mScores),gregexpr("_[0-9]*_",rownames(mScores)))))
ends=gsub('_','',unlist(regmatches(rownames(mScores),gregexpr("_[0-9]*$",rownames(mScores)))))
SNPpos=(as.numeric(starts)+as.numeric(ends)-1)/2
mScores2=cbind(mScores,chr=chromos,SNPpos=SNPpos)
mScores=mScores2
save(mScores,file='/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-06-22/motifScores.acrossTFs.Above0MaxLogodds.withAnno.rda')

###### Now, the motif part can be directly loaded, no need to rerun the above code! #######
#==========================================================================================

