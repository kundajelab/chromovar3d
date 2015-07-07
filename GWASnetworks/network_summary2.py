from optparse import OptionParser
import os
import pickle
import re
import pdb
import math

def main():
	parser=OptionParser()
	
	parser.add_option('--p',dest='p',default='5')
	parser.add_option('--infile',dest='infile',default='/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-05-30//gwasOverlaps/overlapQTL_EUR.CD_pruned_rsq_0.8_expanded_rsq_0.8/overlapQTL_EUR.CD_pruned_rsq_0.8_expanded_rsq_0.8.bed.LD08_p5.bed.overlapGWAS.network.relevant')
	parser.add_option('--peak_snp_dict',dest='di',default='/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-06-09/QTLs//snp2pvalue.2015-06-09.LocalAndDistal.dict')
	opts,args=parser.parse_args()

	pd=pickle.load( open( opts.di, "rb" ) )
	#pdb.set_trace()
	print 'loaded dictionary'

	out=open(opts.infile+'.net.'+'pval'+opts.p,'w')
	#outtable=open(opts.infile+'.consolidatedTable.txt','w')
	outlabels=open(opts.infile+'.net.labels','w')
	outcolors=open(opts.infile+'.net.colors','w')
	#out.write('node1\tnode2\tLocalDistal\tQTLType\n')
	net={}
	labels={}
	colors={} #reg elt, gene or gwas
	for line in open(opts.infile,'r').readlines():
		items=line.strip().split('\t')
		gwasp=float(items[18])
		if gwasp>(float(math.pow(10,(-float(opts.p))))):
			continue
		#print items
		#the qtl snp
		snp=items[3]
		if snp not in labels.keys():
			labels[snp]=snp 
		#the gwas snp
		gwas=items[17]
		if str('GWAS:'+gwas) not in labels.keys():
			labels['GWAS:'+gwas]=gwas
			colors['GWAS:'+gwas]='GWAS'
		#peak names
		local_pname=items[6]
		distal_pname=items[7]
		if items[10] not in colors.keys():
			colors[items[10]]='enh'
		if items[11] not in colors.keys():
			colors[items[11]]='enh'	
		if items[12] not in colors.keys():
			colors[items[12]]='tss'
		if items[13] not in colors.keys():
			colors[items[13]]='tss'	
		#gene assignments
		local_peak=items[12] #gene
		distal_peak=items[13]
		if local_peak not in labels.keys():
			labels[local_peak]=local_peak
		if distal_peak not in labels.keys():
			labels[distal_peak]=distal_peak
		if local_peak=='NA':
		   local_peak=items[10] #this is a reg elt
		if distal_peak=='NA':
		   distal_peak=items[11] #also a reg elt
		#if reg element, then no label
		if local_peak not in labels.keys():
			labels[local_peak]='.'
		if distal_peak not in labels.keys():
			labels[distal_peak]='.'
		#also compute midpoint of distal (we'll use to show distance to qtl)
		if items[11]!='NA':
			distal_coords=re.sub(':','-',items[11]).split('-')
			distal_midpoint=(int(distal_coords[1])+int(distal_coords[2]))/2
			qtl_coords=re.sub(':','-',items[4]).split('-')[2]
		#type of QTL (for coloring edges)
		qtl_type_local='chromatin'
		if 'RNA' in items[6]:
			qtl_type_local='RNA'
		qtl_type_distal='chromatin'
		if 'RNA' in items[7]:
			qtl_type_distal='RNA'
		#local pvalue
		local_p=pd['local'][snp][local_pname]

		

		#make network
		if local_peak not in net.keys():
			net[local_peak]={}
			net[local_peak]['snps']={}
			net[local_peak]['snps'][qtl_type_local]=set()
			net[local_peak]['gwas']=set()

			net[local_peak]['bestSNP']={}
			net[local_peak]['bestSNP'][qtl_type_local]={}
			net[local_peak]['bestSNP'][qtl_type_local]['snp']=snp
			net[local_peak]['bestSNP'][qtl_type_local]['p']=local_p
			net[local_peak]['peak_name']=set() 

		if qtl_type_local not in net[local_peak]['snps'].keys():
			net[local_peak]['snps'][qtl_type_local]=set()
			net[local_peak]['bestSNP'][qtl_type_local]={}
			net[local_peak]['bestSNP'][qtl_type_local]['snp']=snp
			net[local_peak]['bestSNP'][qtl_type_local]['p']=local_p
		#print items[5]
		#print net[local_peak]
		if local_p<net[local_peak]['bestSNP'][qtl_type_local]['p']:
			net[local_peak]['bestSNP'][qtl_type_local]['snp']=snp 
			net[local_peak]['bestSNP'][qtl_type_local]['p']=local_p 
		net[local_peak]['peak_name'].add(local_pname)

		if items[5]=='Distal':
			#pvalue for the snp-peak
			#distal_p=pd['distal'][snp][distal_peak]
			if distal_peak not in net[local_peak].keys(): 
				net[local_peak][distal_peak]={}
				net[local_peak][distal_peak]['snps']={}
				net[local_peak][distal_peak]['snps'][qtl_type_distal]=set() #these will be the snps for this local/distal pair
				net[local_peak][distal_peak]['gwas']=set()
				net[local_peak][distal_peak]['dist']={}
				
				#net[local_peak][distal_peak]['bestSNP']={}
				#net[local_peak][distal_peak]['bestSNP'][qtl_type_distal]={}
				#net[local_peak][distal_peak]['bestSNP'][qtl_type_distal]['snp']=snp
				#net[local_peak][distal_peak]['bestSNP'][qtl_type_distal]['p']=pd['distal'][snp][distal_pname]

			if qtl_type_distal not in net[local_peak][distal_peak]['snps'].keys():
				net[local_peak][distal_peak]['snps'][qtl_type_distal]=set()
			net[local_peak][distal_peak]['snps'][qtl_type_distal].add(snp) #add snp as distal qtl
			net[local_peak][distal_peak]['gwas'].add(gwas)
			net[local_peak][distal_peak]['dist'][snp]=abs(int(distal_midpoint)-int(qtl_coords))/1000
		
			
			
		net[local_peak]['snps'][qtl_type_local].add(snp) #add the snp as a local QTL
		net[local_peak]['gwas'].add(gwas)

	#print net


	#write network
	for local_peak in net.keys():
		gwases=set()
		distals=set(net[local_peak].keys()).difference(set(['gwas','snps','bestSNP','peak_name']))
		snps_added=set()
		for distal_peak in distals:
			qtl_types=net[local_peak][distal_peak]['snps'].keys()
			for qtl_type in qtl_types:
				#get a snp for the distal
				########TODO: more principled way to pick the snp
				#print 'Distal choosing from: '+str(len(list(net[local_peak][distal_peak]['snps'][qtl_type])))
				#print list(net[local_peak][distal_peak]['snps'][qtl_type])
				p=float(1)
				distalsnp=net[local_peak][distal_peak]['snps'][qtl_type]
				for candidate in list(net[local_peak][distal_peak]['snps'][qtl_type]):
					#print candidate
					for local_peak_name in net[local_peak]['peak_name']:
						if local_peak_name in pd['local'][candidate].keys():
							print pd['local'][candidate]
							cur_p=float(pd['local'][candidate][local_peak_name])
							if cur_p<p:
								p=cur_p
								distalsnp=candidate
				#print 'Picked '+str(distalsnp)
				##as it was
				'''
				distalsnp=list(net[local_peak][distal_peak]['snps'][qtl_type])[0]
				print 'Distal choosing from: '+str(len(list(net[local_peak][distal_peak]['snps'][qtl_type])))
				print list(net[local_peak][distal_peak]['snps'][qtl_type])
				'''
				#distalsnp=net[local_peak][distal_peak][qtl_type_distal]['bestSNP']['snp']
				snps_added.add(distalsnp)
				out.write(distalsnp+'\t'+distal_peak+'\t'+'Distal'+'\t'+qtl_type+'\t'+str(net[local_peak][distal_peak]['dist'][distalsnp])+' kb'+'\n')
				#remember the gwases for this one
				for gwas in net[local_peak][distal_peak]['gwas']:
					out.write('GWAS:'+gwas+'\t'+distalsnp+'\t'+'LD'+'\t'+'NA\t.'+'\n')
				#write the local now connected to the same snp as the distal
				#- first figure out the type QTL we had for the local
				for possible_local_qtl_type in net[local_peak]['snps'].keys():
					if distalsnp in net[local_peak]['snps'][possible_local_qtl_type]:
						out.write(distalsnp+'\t'+local_peak+'\t'+'Local'+'\t'+possible_local_qtl_type+'\t.'+'\n')
						#outtable.write(gwas+'\t'+distalsnp+'\t'+local_peak+'\t'+'Local'+'\t'+possible_local_qtl_type+'\t.'+'\n')
		
		for gwas in net[local_peak]['gwas']:
			#connect it it any leftover local snp
			if len(distals)==0: #there were no distal QTLs, so no chances to add a QTL snp
				##########TODO: more principled way to pick the snp
				for possible_local_qtl_type in net[local_peak]['snps'].keys():
					snp=net[local_peak]['bestSNP'][possible_local_qtl_type]['snp']
					#print 'Local. Picked best SNP' 
					#list(net[local_peak]['snps'][possible_local_qtl_type])[0]
					out.write('GWAS:'+gwas+'\t'+snp+'\t'+'LD'+'\t'+'NA\t.'+'\n')
					if snp not in snps_added:
						out.write(snp+'\t'+local_peak+'\t'+'Local'+'\t'+possible_local_qtl_type+'\t.'+'\n')
						#outtable.write(gwas+'\t'+snp+'\t'+local_peak+'\t'+'Local'+'\t'+possible_local_qtl_type+'\t.'+'\n')
	#write node labels
	outlabels.write('node\tlabel\n')
	outcolors.write('node\tcolor\n')
	for node in labels.keys():
		outlabels.write(node+'\t'+labels[node]+'\n')
	for node in colors.keys():
		outcolors.write(node+'\t'+colors[node]+'\n')

	#print opts.infile+'.small'

main()