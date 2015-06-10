from optparse import OptionParser
import os
import math,re
from time import gmtime, strftime
import gzip,sys
'''
Make vcf for imputed individuals.
'''

def main():
    parser=OptionParser()
    parser.add_option('--mfile',dest='mfile',default='/srv/gsfs0/projects/kundaje/commonRepository/encode/data/byDataType/bindingSites/encodeMotifs/dec2013/motifs/motifs.txt')
    parser.add_option('--out',dest='out',default='/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/data/motifPWM/pouya.PWM.motif.Homer')
    parser.add_option('--out_match',dest='out_match',default='/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/data/motifPWM/matchesPerTF/')
    parser.add_option('--thresh',dest='thresh',default='4')
    opts,args=parser.parse_args()
    
    os.system('mkdir '+opts.out_match)
    out=open(opts.out,'w')
    for line in open(opts.mfile,'r').readlines():
        if line[0:1]=='>':
            try:
                out_cur.close()
            except:
                print 'No cur file open'
            out_cur=open(opts.out+re.sub('>','.Motif.',line.strip().split()[0]),'w')
            out_cur_forscan=open(opts.out+re.sub('>','.Motif.',line.strip()).split()[0]+'scanThresh'+opts.thresh,'w')
            #This is a motif name
            out.write('>Consensus'+'\t'+re.sub(' ','_',line.strip())+'\t'+'0'+'\n')
            out_cur.write('>Consensus'+'\t'+re.sub(' ','_',line.strip())+'\t'+'0'+'\n')
            out_cur_forscan.write('>Consensus'+'\t'+re.sub('\t','',re.sub(' ','_',line.strip()))+'\t'+opts.thresh+'\n')
            #Figure out the mptif name and get its matches in the genome
            motif_name=re.sub('>','',re.sub(' ','_',line.strip()).split()[0])
            print motif_name
            if 'known1' in motif_name:
                #if 'SPI1' not in motif_name:
                #    continue
                cmd_motifmatch='zcat '+opts.motif_match+' | grep '+motif_name+' | cut -d \' \' -f2-4 | gzip > '+opts.out_match+os.path.basename(opts.out)+re.sub('>','.Motif.',line.strip().split()[0])+'.MATCH.bed.gz'
                print cmd_motifmatch
                #qsub_a_command(cmd_motifmatch,opts.out_match+os.path.basename(opts.out)+re.sub('>','.Motif.',line.strip().split()[0])+'.MATCH.bed.gz'+'_script.sh','qqqq','3G')
        else:
            #This is from PWM. We add 0.001 and we write it.
            items=line.strip().split()
            psc=0.001
            out.write(str(psc+float(items[1]))+'\t'+str(psc+float(items[2]))+'\t'+str(psc+float(items[3]))+'\t'+str(psc+float(items[4]))+'\n')
            out_cur.write(str(psc+float(items[1]))+'\t'+str(psc+float(items[2]))+'\t'+str(psc+float(items[3]))+'\t'+str(psc+float(items[4]))+'\n')
            out_cur_forscan.write(str(psc+float(items[1]))+'\t'+str(psc+float(items[2]))+'\t'+str(psc+float(items[3]))+'\t'+str(psc+float(items[4]))+'\n')

def qsub_a_command(cmd,shell_script_name,split_string=',',memory_number='20G'):
    f=open(shell_script_name,'w')
    cmds=cmd.split(split_string)
    for i in range(len(cmds)):
        f.write(cmds[i]+'\n')
    f.close()
    os.system('chmod 711 '+shell_script_name)
    #Qsub the script                                                                                                                                                       
    os.system("qsub -l mem_free="+memory_number+" -l h_vmem="+memory_number+" -l h_rt=20:00:00 -o "+shell_script_name+'.o'+' -e '+shell_script_name+'.e'+' '+shell_script_name)



                


main()
