OUTNAME='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-05-30'
DISTAL='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-05-30/CrossTFanalysis/DistalAnalysis.DistalPeakIs.gz.CompiledDataAcrossTFs'
LOCAL='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-05-30/CrossTFanalysis/LocalAnalysis.motifMAF05..gz.annotated.gz.CompiledDataAcrossTFs'

distal=read.table(DISTAL,header=TRUE)
local=read.table(LOCAL,header=TRUE)

distal=distal[,c('snp',
	'affected.gene',
	'snp.pos',
	'localpeak',
	'spearmanCorr_MotifScoreVSpeak',
	'maxObsMotif',
	'motifName',
	'genomeWide_match.p.sp',
	'affected.hmark','crossTF_pBH.spearman')]
colnames(distal)=gsub('affected.gene','affected.peakORgene',colnames(distal))
local=local[,c('SNP',
	'affected.gene',
	'snp.position',
	'localpeak',
	'spearmanCorr_MotifScoreVSpeak',
	'maxObsMotif',
	'motifName',
	'genomeWide_match.p.sp',
	'affected.hmark','crossTF_pBH.spearman')]
colnames(local)=gsub('affected.gene','affected.peakORgene',
	gsub('SNP','snp',
		gsub('snp.position','snp.pos',
			colnames(local))))
local=cbind(local,LocalOrDistal='Local')
distal=cbind(distal,LocalOrDistal='Distal')

motif=rbind(local,distal)

write.table(motif,
	file=paste(OUTNAME,'/MotifAnalysisTable.txt',sep=''),
	quote=F,sep='\t',row.names=F,col.names=T)