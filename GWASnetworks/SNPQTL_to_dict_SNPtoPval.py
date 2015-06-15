from optparse import OptionParser
import os
import subprocess
import re
import pickle

def main():
	parser=OptionParser()
	
	parser.add_option('--snpqtl_local',dest='snpqtl_local',default='/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-06-09/QTLs//pvalueSortedSNPspeakPeak.Local')
	parser.add_option('--snpqtl_distal',dest='snpqtl_distal',default='/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-06-09/QTLs//pvalueSortedSNPspeakPeak.Distal')
	parser.add_option('--out',dest='out',default='/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-06-09/QTLs//snp2pvalue.2015-06-09.LocalAndDistal.dict')
	opts,args=parser.parse_args()

	local=get_dict_from_file(opts.snpqtl_local)
	distal=get_dict_from_file(opts.snpqtl_distal)
	d={}
	d['local']=local
	d['distal']=distal 
	
	pickle.dump( d, open( opts.out, "wb" ) )
	print opts.out

def get_dict_from_file(f):
	d={}
	c=0
	for line in open(f,'r').readlines():
		print c
		c+=1
		pos,snp,pval,peak=line.strip().split('\t')
		if snp not in d.keys():
			d[snp]={}
		if peak not in d[snp].keys():
			d[snp][peak]=pval
	return d
main()