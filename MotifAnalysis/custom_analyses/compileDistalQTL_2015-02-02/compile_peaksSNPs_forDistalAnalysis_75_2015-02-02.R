
DISTAL_QTLS='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Distal/DistalQTLs.2015-02-02'

#Add peak annotations to the distal QTL file
#============================================
#one line for each peak combo
hmarks=c('H3K4ME1','H3K4ME3','H3K27AC','dhs','RNA')
first=TRUE
for (hmark in hmarks){
    qtl.f=paste('/srv/gs1/projects/snyder/jzaugg/histoneQTL/hQTL_analysis/results/20140622_newPeaks/allQTLs_states.prox.frag.FDR10.allChr.2000kb.',hmark,'.HiC_TH0.4.LDcorrected_0.2.states.rda',sep='')
    load(qtl.f)
    allQTLs=allQTLs_states
    allQTLs$distalQTL.peak.id=paste(hmark,'_',allQTLs$distalQTL.peak.id,sep='')
    if (first==FALSE) {
       qtl=rbind(qtl,allQTLs)
    }
    if (first==TRUE){
       qtl=allQTLs
       first=FALSE
    }
}
#go through each line of the qtl file and make a new dataset
#we are making a separate line for each local QTL that can influence a distalQTL
qtl_singleRow=qtl[1,]
qtl_singleRow=qtl_singleRow[-1,]
for (i in c(1:nrow(qtl))){
    print(i)
    locals=strsplit(qtl[i,'localQTL.peak.id'],',')[[1]]
    for (local_peak in locals){
    	cur_row=qtl[i,]
	    cur_row$localQTL.peak.id=local_peak
    	qtl_singleRow=rbind(qtl_singleRow,cur_row)
    }
}
COLS_DISCARD=c('distalQTL_enh.state','localQTL_enh.state',
  'distalQTL_tss.state','localQTL_tss.state',
  'distalQTL_enh.tss.state','localQTL_enh.tss.state')
qtl_singleRow_final=qtl_singleRow[,-which(colnames(qtl_singleRow) %in% COLS_DISCARD)]
#write down this dataset of all distalQTLs combined
write.table(qtl_singleRow_final,file=DISTAL_QTLS,sep='\t',quote=F,row.names=F,col.names=T)

#Overlap RNA with any peak1, and add it as a testable peak1.
#===========================================================
require(rtracklayer)
require(GenomicRanges)
rnaanno=read.table('/srv/gs1/projects/snyder/jzaugg/histoneQTL/peakAnnotation/peakLocations.ALL.RNA.txt',header=T)
rna.range=GRanges(seqnames =Rle(rnaanno[,'chr']), ranges =IRanges(rnaanno[,'start'], end = rnaanno[,'end'])) #rna
newqtl=qtl_singleRow

newqtl_brief=newqtl[,c('chr','distalQTL.peak.id','localQTL.peak.id')]
newqtl_with_overlappingRNA=newqtl_brief
unique_peak1s=as.character(unique(newqtl_brief$distalQTL.peak.id))
pcounter=0
for (p1 in unique_peak1s){
    pcounter=pcounter+1
    print(pcounter)
    p1.rows=newqtl[which(as.character(newqtl$distalQTL.peak.id)==p1),]
    p1.range=GRanges(seqnames =Rle(p1.rows[1,'chr']), ranges =IRanges(p1.rows[1,'peak.start'], end = p1.rows[1,'peak.end']))
    matching.rnas=as.data.frame(findOverlaps(query=rna.range,subject=p1.range,select='all'))[,'queryHits']
    if (length(matching.rnas)>0){
       for (rna.item in c(1:length(matching.rnas))){
       	   rna.names=paste('RNA_',rnaanno[as.numeric(as.character(matching.rnas[rna.item])),'id'],sep='')
       	   combinedRows=p1.rows[,c('chr','distalQTL.peak.id','localQTL.peak.id')]
       	   combinedRows[,'distalQTL.peak.id']=rna.names
	   print(combinedRows)
       	   newqtl_with_overlappingRNA=rbind(newqtl_with_overlappingRNA,combinedRows)
       }
    }
}

write.table(newqtl_with_overlappingRNA,file='/home/oursu/rna')
#remove duplicated rows too
newqtl_with_overlappingRNA.dupli=which(duplicated(newqtl_with_overlappingRNA))
if (length(newqtl_with_overlappingRNA.dupli)>0){
   newqtl_with_overlappingRNA=newqtl_with_overlappingRNA[-newqtl_with_overlappingRNA.dupli,]
}

newqtl2=newqtl[,c('chr','distalQTL.peak.id','localQTL.peak.id')]
newqtl2.uniq=newqtl2[-which(duplicated(newqtl2)),]
write.table(newqtl2.uniq,file=paste(DISTAL_QTLS,'.ChrGeneDistal.txt',sep=''),
  sep='\t',quote=F,row.names=F,col.names=T)
write.table(newqtl_with_overlappingRNA,file=paste(DISTAL_QTLS,'.ChrGeneDistal.PlusOverlappedRNA.txt',sep=''),
  sep='\t',quote=F,row.names=F,col.names=T)
