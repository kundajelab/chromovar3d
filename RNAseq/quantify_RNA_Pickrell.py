from optparse import OptionParser
import os,sys
import math
from time import gmtime, strftime

'''
Author:Oana Ursu
'''

def main():
    parser=OptionParser()
    
    parser.add_option('--out',dest='out',help='Out',default='/srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/results/Pickrell/SAILFISH/')
    our_study='/srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.ALL_for_peakCalling.individuals_to_keep'
    pickrell='/srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/data/list_lanes_pickrell_2010_nature'
    opts,args=parser.parse_args()
    
    our_people=set()
    for line in open(our_study,'r').readlines():
        person=line.strip().split('\t')[0]
        our_people.add(person)
    of_interest=our_people
    gender={}
    for line in open(our_study,'r').readlines():
        our_people.add(line.strip().split('\t')[0])
        items=line.strip().split('\t')
        person=items[0]
        person_gender=items[5]
        if person not in gender.keys():
            gender[person]=person_gender
      
    pickrell_dir='/srv/gsfs0/projects/kundaje/users/oursu/histoneQTLproject/RNAseq/results/Pickrell/'

    person_to_fq={}
    dirs_in_pickrell=os.listdir(pickrell_dir)
    for f in dirs_in_pickrell:
        if 'NA' in f:
            person=f
            if person not in of_interest:
                continue
            if person not in person_to_fq.keys():
                person_to_fq[person]=set()
            curr_fq=pickrell_dir+person+'/'+person+'_argonne.fastq.gz'
            if os.path.isfile(curr_fq):
                person_to_fq[person].add(curr_fq)
            curr2=pickrell_dir+person+'/'+person+'_2_argonne.fastq.gz'
            if os.path.isfile(curr2):
                person_to_fq[person].add(curr2)
    for i in person_to_fq.keys():
        print person_to_fq[i]
    
    
    for person in person_to_fq.keys():
        out_curr=opts.out+person+'.sailfish'
        #And run sailfish
        cur_gender=gender[person]
        sailfish_idx='/srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/data/Transcriptome/gencode.v19.annotation.PC.lincRNA.gtf.splicedExon.N'+cur_gender+'.fa.dedup.fa_IDX_sailfish'
        cmd_module='module load sailfish/0.6.3'
        library_type='"T=SE:S=U"' #T=PE:O=><:S=SA
        fastqs=list(person_to_fq[person])
        print fastqs
        cmds=[]
        cmds.append('#!/usr/bin/env bash')
        cmds.append(cmd_module)
        sailfish_cmd='sailfish quant -i '+sailfish_idx+' -l '+library_type+' -r <(gunzip -c '+' '.join(fastqs)+') -o '+out_curr+' -f'
        cmds.append(sailfish_cmd)
        cmds.append('cd '+out_curr)
        cmds.append('module load java/latest')
        gtf='/srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/data/GENCODE_v19_2014-06-03/gencode.v19.annotation.PC.lincRNA.gtf'
        cmds.append('/srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/src/TranscriptsToGenes.sh --gtf-file '+gtf+' --exp-file '+out_curr+'/quant.sf'+' --res-file '+person+'quant.gene_level.sf')
        cmds.append('mv '+out_curr+'/quant.sf'+' '+out_curr+'/'+person+'quant.sf')
        print '\n'.join(cmds)
        qsub_a_command('qqqq'.join(cmds),out_curr+'_script.sh','qqqq','10G')

def qsub_a_command(cmd,shell_script_name,split_string=',',memory_number='20G'):
    f=open(shell_script_name,'w')
    cmds=cmd.split(split_string)
    for i in range(len(cmds)):
        f.write(cmds[i]+'\n')
    f.close()
    os.system('chmod 711 '+shell_script_name)
    os.system("qsub -l mem_free="+memory_number+" -l h_vmem="+memory_number+" -l h_rt=20:00:00 -o "+shell_script_name+'.o'+' -e '+shell_script_name+'.e'+' '+shell_script_name)
    


main()
