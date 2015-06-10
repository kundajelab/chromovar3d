
require(rtracklayer)
require(GenomicRanges)

i.f='/srv/gs1/projects/snyder/jzaugg/histoneQTL/HiC/InteractionPipeline/GM12878_Fabian/SignificantInteractions/HiCcorrelations_0.2_2MB_chrCHROM.txt'
chromo=1
interactions=read.table(gsub('CHROM',as.character(chromo),i.f),header=T,nrow=100)
countPos=table(c(interactions$pos_i,interactions$pos_j))
iCount=data.frame(count=countPos,fragName=names(countPos),chr=chromo)
for (chromo_idx in c(2:22)){
    interactions_cur=read.table(gsub('CHROM',as.character(chromo_idx),i.f),header=T,nrow=100)
    countPos=table(c(interactions_cur$pos_i,interactions_cur$pos_j))
    iCount=rbind(iCount,data.frame(count=countPos,fragName=paste('chr',chromo_idx,'_',names(countPos),sep=''),chr=chromo_idx))	
}

#Read in distal
qtl.f='/srv/gs1/projects/snyder/jzaugg/histoneQTL/hQTL_analysis/results/20140622_newPeaks/DistQTL_list/HiC_local_QTL_dist_QTL.all.txt'
qtl=read.table(qtl.f,header=T,sep='\t')
qtl_peak1= GRanges(seqnames =Rle(paste('chr',qtl[,'chr'],sep='')), ranges =IRanges(qtl[,'peak1.start'], end = qtl[,'peak1.end']))
qtl_peak2= GRanges(seqnames =Rle(paste('chr',qtl[,'chr'],sep='')), ranges =IRanges(qtl[,'peak2.start'], end = qtl[,'peak2.end']))

#Now, make a ranges for the fragments as well
fragments.f='/srv/gsfs0/projects/kundaje/users/oursu/linkingChromatin/old/data/2014-05-11/HiC_fragments.bed'
frags=read.table(fragments.f)
frags=cbind(frags,mid=paste(frags[,1],'_',(frags[,2]+frags[,3])/2,sep=''))
frag_range=GRanges(seqnames =Rle(frags[,1]), ranges =IRanges(frags[,2], end = frags[,2]))

frags_in_peaks=as.data.frame(findOverlaps(query=qtl_peak1,subject=frag_range,select='all'))
colnames(frags_in_peaks)=c('peak1','frag')
peak2frag=cbind(peak1=rownames(as.data.frame(qtl_peak1))[frags_in_peaks[,'peak1']],frag=rownames(as.data.frame(frag_range))[frags_in_peaks[,'frag']])