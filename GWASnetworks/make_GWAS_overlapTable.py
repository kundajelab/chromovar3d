from optparse import OptionParser
import os
import pickle
import re
import pdb

def main():
	parser=OptionParser()
	
	parser.add_option('--infile',dest='infile',default='/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-06-13/gwasOverlaps/overlapQTL_EUR.CD_pruned_rsq_0.8_expanded_rsq_0.8/overlapQTL_EUR.CD_pruned_rsq_0.8_expanded_rsq_0.8.bed.LD08_p5.bed.overlapGWAS.bed.net.relevant')
	parser.add_option('--peak_snp_dict',dest='di',default='/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-06-13/QTLs//snp2pvalue.2015-06-09.LocalAndDistal.dict')
	opts,args=parser.parse_args()

	pd=pickle.load( open( opts.di, "rb" ) )
	print 'loaded dictionary'
	#pdb.set_trace()
	out=open(opts.infile+'.table.txt','w')
	#net=open(opts.infile+'.network.txt','w')

	DELIMITER=';'
	d={}
	for line in open(opts.infile,'r').readlines():
		items=line.strip().split('\t')
		#the data
		chrLDSNP,startLDSNP,endLDSNP,QTLSNP,QTLSNPcoords,QTLType,LocalPeak,DistalPeak,R2,LDSNP,RegEltLocal,RegEltDistal,TSSLocal,TSSDistal,GWASchr,GWASstart,GWASend,GWAStagSNPID,GWASassociationPval=items
		gwas_coords=GWASchr+':'+GWASstart+'-'+GWASend
		local_qtl_type='chromatin'
		distal_qtl_type='chromatin'
		if 'RNA' in LocalPeak:
			local_qtl_type='RNA'
		if 'RNA' in DistalPeak:
			distal_qtl_type='RNA'
		local_qtl_type='NA'
		distal_qtl_type='NA'
		#print items
		l_d_combo='LocalElt:'+RegEltLocal+',DistalElt:'+RegEltDistal#+',LocalGene:'+TSSLocal+',DistalGene:'+TSSDistal
		if l_d_combo not in d.keys():
			d[l_d_combo]={}
			#gwas data
			d[l_d_combo]['GWAS']={}
			d[l_d_combo]['RegEltLocal']=RegEltLocal
			d[l_d_combo]['RegEltDistal']=RegEltDistal
			d[l_d_combo]['TSSLocal']=TSSLocal
			d[l_d_combo]['TSSDistal']=TSSDistal
			#Local data
			d[l_d_combo]['Local']={}
			d[l_d_combo]['Local']['Peaks']=set()
			d[l_d_combo]['Local']['Genes']=set()
			#Distal data
			d[l_d_combo]['Distal']={}
			d[l_d_combo]['Distal']['Peaks']=set()
			d[l_d_combo]['Distal']['Genes']=set()

		#gwas data
		if GWAStagSNPID not in d[l_d_combo]['GWAS'].keys():
			d[l_d_combo]['GWAS'][GWAStagSNPID]={}
			#QTL SNP
			d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']={}
			d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local']={}
			d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal']={}

		#update gwas data
		d[l_d_combo]['GWAS'][GWAStagSNPID]['coords']=gwas_coords
		d[l_d_combo]['GWAS'][GWAStagSNPID]['p']=GWASassociationPval
		#d[l_d_combo]['GWAS'][GWAStagSNPID]['tagSNP']=GWAStagSNPID

		#update local data
		d[l_d_combo]['Local']['Peaks'].add(LocalPeak)
		d[l_d_combo]['Local']['Genes'].add(TSSLocal)

		#update distal data
		d[l_d_combo]['Distal']['Peaks'].add(DistalPeak)
		d[l_d_combo]['Distal']['Genes'].add(TSSDistal)

		#QTL SNP
		#Local
		if local_qtl_type not in d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'].keys():
			d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type]={}
		if QTLSNP not in d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type].keys():
			d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type][QTLSNP]={}
			d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type][QTLSNP]['LD']={}
			d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type][QTLSNP]['p']=1
			d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type][QTLSNP]['peak']='TBA'
		if gwas_coords not in d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type][QTLSNP]['LD'].keys():
			d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type][QTLSNP]['LD'][gwas_coords]=R2
		p=float(pd['local'][QTLSNP][LocalPeak])
		if p<d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type][QTLSNP]['p']:
			d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type][QTLSNP]['p']=p
			d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type][QTLSNP]['peak']=LocalPeak
		#Distal
		#only update if it's a distal QTL
		if QTLType=='Distal':
			if distal_qtl_type not in d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'].keys():
				d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type]={}
			if QTLSNP not in d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type].keys():
				d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type][QTLSNP]={}
				d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type][QTLSNP]['LD']={}
				d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type][QTLSNP]['p']=1
				d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type][QTLSNP]['peak']='TBA'
			if gwas_coords not in d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type][QTLSNP]['LD'].keys():
				d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type][QTLSNP]['LD'][gwas_coords]=R2
			p=float(pd['distal'][QTLSNP][DistalPeak])
			#print p
			if p<d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type][QTLSNP]['p']:
				d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type][QTLSNP]['p']=p
				#print 'putting in peak'
				#print DistalPeak
				d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type][QTLSNP]['peak']=DistalPeak

	#Now write down this network
	net_lines=set()
	for l_d_combo in d.keys(): #for each pair
		for GWAStagSNPID in d[l_d_combo]['GWAS'].keys(): #we get 1 line per gwas tag snp
			
			#write network====================
			'''
			for qtl_type in ['chromatin','RNA']:
				#- GWAS tag SNP---->QTL SNP with best pvalue
				best_local=get_best_qtl(d[l_d_combo]['GWAS'][GWAStagSNPID],'Local',qtl_type)
				if best_local=='TBA':
					continue
				net_lines.add(GWAStagSNPID+'\tLD\t'+best_local+'.'+'NA')
				#- QTL SNP ----> local
				#try to connect to the regElt first
				if d[l_d_combo]['RegEltLocal']!='NA':
					net_lines.add(best_local+'\tLocal\t'+d[l_d_combo]['RegEltLocal']+'\t'+'.'+'\t'+qtl_type)
				#then try tss
				if d[_l_combo]['TSSLocal']!='NA':
					net_lines.add(best_local+'\tLocal\t'+d[l_d_combo]['TSSLocal']+'\t'+'.'+'\t'+qtl_type)
			'''

			#write table======================
			outtext=GWAStagSNPID+'\t'
			outtext+=d[l_d_combo]['GWAS'][GWAStagSNPID]['p']+'\t'
			outtext+=l_d_combo+'\t'
			outtext+=','.join(list(d[l_d_combo]['Local']['Peaks']))+'\t'
			outtext+=','.join(list(d[l_d_combo]['Local']['Genes']))+'\t'
			outtext+=','.join(list(d[l_d_combo]['Distal']['Peaks']))+'\t'
			outtext+=','.join(list(d[l_d_combo]['Distal']['Genes']))+'\t'
			#Local SNP
			localsnp=''
			first=True
			for local_qtl_type in d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'].keys():
				for QTLSNP in d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type].keys():
					#for i in range(len(d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type][QTLSNP]['p'])):
					if True:
						if not first:
							localsnp+='|'
						localsnp+=QTLSNP+DELIMITER
						localsnp+=d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type][QTLSNP]['peak']+DELIMITER
						localsnp+=str(d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type][QTLSNP]['p'])+DELIMITER
						first2=True
						for gwas_coords in d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type][QTLSNP]['LD'].keys():
							if not first2:
								localsnp+=','
							first2=False
							localsnp+='LD='+gwas_coords+' '+d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][local_qtl_type][QTLSNP]['LD'][gwas_coords]
							#print 'LD='+gwas_coords+' '+d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Local'][QTLSNP]['LD'][gwas_coords]
						localsnp+=DELIMITER
						first=False
			outtext+=localsnp+'\t' 
			#Distal SNP
			distalsnp=''
			if len(d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'].keys())==0:
				#then there are no distals
				outtext+='NA' #for the distalSNP
				print outtext
				out.write(outtext+'\n')
				continue
			first=True
			for distal_qtl_type in d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'].keys():
				for QTLSNP in d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type].keys():
					#for i in range(len(d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type][QTLSNP]['p'])):
					if True:
						if not first:
							distalsnp+='|'
						distalsnp+=QTLSNP+DELIMITER
						distalsnp+=d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type][QTLSNP]['peak']+DELIMITER
						distalsnp+=str(d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type][QTLSNP]['p'])+DELIMITER
						first2=True
						for gwas_coords in d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type][QTLSNP]['LD'].keys():
							if not first2:
								localsnp+=','
							first2=False
							distalsnp+='LD='+gwas_coords+' '+d[l_d_combo]['GWAS'][GWAStagSNPID]['QTL_SNP']['Distal'][distal_qtl_type][QTLSNP]['LD'][gwas_coords]
						distalsnp+=DELIMITER
						first=False
				outtext+=distalsnp 
			out.write(outtext+'\n')
			print outtext
			#print "----------"
	#print d	

def get_best_qtl(d,localdistal,qtl_type):
	snp='TBA'
	p=1
	if qtl_type not in d['QTL_SNP'][localdistal].keys():
		return snp
	for QTLSNP in d['QTL_SNP'][localdistal][qtl_type].keys():
		cur_p=d['QTL_SNP'][localdistal][qtl_type][QTLSNP]['p']
		if cur_p<p:
			snp=QTLSNP 
	return snp










main()