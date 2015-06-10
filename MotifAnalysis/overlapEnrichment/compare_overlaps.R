TF.f='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/OverlapEnrichment/TFBS_overlap_QTLpeaks0kb/TFBS_overlap_H3K27ACQTLpeaks0kb___.overlapMatrix.gz'
motif.f='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/OverlapEnrichment/Motif_overlap_QTLpeaks0kb/Motif_overlap_H3K27ACQTLpeaks0kb___.overlapMatrix.gz'
corMotif.f='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/OverlapEnrichment/MotifCorrelatedLocal_overlap_QTLpeaks0kb/MotifCorrelatedLocal_overlap_H3K27ACQTLpeaks0kb___.overlapMatrix.gz'

fs=c(TF.f,motif.f,corMotif.f)

peaks=list()
tf='CTCF'
for (f in fs){
	data=read.table(f,header=TRUE)
	tfcol=which(grepl(tf,colnames(data)))
	peaks[[basename(f)]]=data[which(as.numeric(as.character(data[,tfcol]))>0),'peak']
}

items=names(peaks)
h=data.frame(array(0,dim=c(length(items),length(items))))
rownames(h)=colnames(h)=items
for (i in items){
	for (j in items){
		h[i,j]=length(intersect(peaks[[i]],peaks[[j]]))/length(peaks[[i]])
	}
}
rownames(h)=colnames(h)=gsub('overlapMatrix.gz','',items)

pdf('/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/OverlapEnrichment/EnrichmentAcrossDatasets/CTCF.pdf')
require(pheatmap)
pheatmap(h,cluster_cols=F,cluster_rows=F,display_numbers=T)
dev.off()