from optparse import OptionParser
import os
import subprocess
import re

def main():
	parser=OptionParser()
	
	parser.add_option('--gtf',dest='gtf',default='/srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/data/GENCODE_v19_2014-06-03/gencode.v19.annotation.PC.lincRNA.gtf')
	parser.add_option('--bed',dest='bed',default='/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-06-09/Annotations/gencode.v19.annotation.PC.lincRNA.TSS.bed')
	opts,args=parser.parse_args()

	out=open(opts.bed,'w')

	for line in open(opts.gtf,'r').readlines():
		items=line.strip().split('\t')
		element=items[2]
		if element!='gene':
			continue
		chromo,start,end,strand=items[0],items[3],items[4],items[6]
		chromo=re.sub('chr','',chromo)
		geneID=re.sub('"','',re.sub('gene_id "','',items[8].split(';')[0]))
		geneName=re.sub('"','',re.sub('gene_name "','',items[8].split(';')[4]))
		if strand=='+':
			tss=max(0,int(start)-1)
		if strand=='-':
			tss=max(0,int(end)-1)

		if len(geneName.strip())==0:
			geneName='NotAvailableGeneSymbol'+geneID
		out.write(chromo+'\t'+str(tss)+'\t'+str(int(tss)+1)+'\t'+geneName+'\n')

main()