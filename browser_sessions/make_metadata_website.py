from optparse import OptionParser
import os
import subprocess
import re
from os import listdir

def main():
	parser=OptionParser()
	
	parser.add_option('--web',dest='web',default='http://chromovar3d.stanford.edu/histone/results/signal/bigwig/')
	parser.add_option('--bigwig_data',dest='bigwig',default='/srv/gs1/projects/snyder/jzaugg/histoneQTL/chromovar3d/histone/results/signal/bigwig/')
	parser.add_option('--out',dest='out',default='/home/oursu/testmetadata')
	parser.add_option('--DNase',dest='DNase',action='store_true')
	opts,args=parser.parse_args()

	out=open(opts.out,'w')

	data={}

	items=listdir(opts.bigwig)
	for item in items:
		if 'bigwig' in item:
			print item
			item_parts=re.sub('.mergedReplicates.subsampleTo50000000.wiggler.norm5.rawsignal.bigwig','',re.sub('SNYDER_','',item)).split('_')
			person=item_parts[1]
                        hmark=item_parts[2]
			if opts.DNase:
				item_parts=re.sub('.subsampleTo35000000.wiggler.norm5.rawsignal.bigwig','',re.sub('SNYDER_','',item)).split('_')
				person=item_parts[2]
				hmark=item_parts[1]
			out.write(person+'\t'+hmark+'\t'+opts.web+item+'\n')
	
	

main()
