from optparse import OptionParser
import os
import subprocess
import re

def main():
	parser=OptionParser()
	
	parser.add_option('--infile',dest='infile',default='/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-06-09/QTLs/chr10.LocalQTLs.bed')
	parser.add_option('--LDfile',dest='LDfile',default='/srv/gsfs0/projects/kundaje/users/pgreens/data/LD/YRI_geno.txt.gz')
	opts,args=parser.parse_args()

	m='module load tabix/0.2.6'
	os.system(m)
	#snps=''
	outf=open(re.sub('.bed','',opts.infile)+'.expandLD_'+re.sub('.gz','',os.path.basename(opts.LDfile))+'.new.bedbedbed','w')

	c=0
	for line in open(opts.infile,'r').readlines():
		print c
		c+=1
		chromo,start,end,name=line.strip().split('\t')
		rsid=name.split(',')[0]
		chromo=re.sub('chr','',chromo)
		lookup_name='chr'+chromo+':'+start+'-'+end

		#include the SNP itself
		outf.write(line.strip()+'\t'+'chr'+chromo+'\t'+start+'\t'+end+'\t'+rsid+'\t'+'chr'+chromo+'\t'+start+'\t'+end+'\t'+rsid+'\t'+'1'+'\n')

		#get snps in LD
		lookup_cmd='tabix '+opts.LDfile+' '+lookup_name #a single locus
		proc = subprocess.Popen(lookup_cmd, stdout=subprocess.PIPE, shell=True)
		(out, err) = proc.communicate()
		for snp in out.split('\n'):
			if len(snp)==0: #not an entry, just the end of the output
				continue 
			if rsid!=snp.strip().split('\t')[3]:
				continue
			outf.write(line.strip()+'\t'+snp+'\n')
		#snps=snps+' '+lookup_name
	#full_cmd='tabix '+opts.LDfile+' '+snps
	

main()