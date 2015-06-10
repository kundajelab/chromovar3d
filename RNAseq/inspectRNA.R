
rna.f='/srv/gsfs0/projects/kundaje/users/oursu/testAlignments/testFastq.fq_2/quant.sf'

rna.f='/srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/results/SAILFISH/NA18520.sailfish/quant.sf'
rna=read.table(rna.f)
colnames(rna)=c('Transcript','Length','TPM','RPKM','KPKM','EstimatedNumKmers','EstimatedNumReads')
pdf(paste(rna.f,'_.pdf',sep=''))
for (i in c(2:ncol(rna))){
    hist(rna[,i],main=colnames(rna)[i],breaks=1000)
    hist(log(rna[,i]),main=paste('log ',colnames(rna)[i]),breaks=100)
}
dev.off()
