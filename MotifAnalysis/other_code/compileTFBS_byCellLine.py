from optparse import OptionParser
import os
import math,re
from time import gmtime, strftime
import gzip,sys
'''
Compile TFBS for cell line
'''

def main():
    parser=OptionParser()
    parser.add_option('--tfbs_dir',dest='tfbsdir',default='/srv/gsfs0/projects/kundaje/commonRepository/encode/data/byDataType/peaks_spp/mar2012/distinct/idrOptimalBlackListFilt/')
    parser.add_option('--out',dest='out',default='/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/data/TFBS_GM12878.bed.gz')
    opts,args=parser.parse_args()
    
    #Find files with Gm12878, Tfbs and no Pol2
    tfbs_files=os.listdir(opts.tfbsdir)
    tfbs_keep=set()
    for tfbs in tfbs_files:
        if 'Pol2' in tfbs:
            continue
        if 'Tfbs' not in tfbs:
            continue
        if 'Gm12878' not in tfbs:
            continue
        print tfbs
        tfbs_keep.add(opts.tfbsdir+tfbs)
    cmds=[]
    cmds.append('module load bedtools/2.18.0')
    cmd='zcat '+' '.join(list(tfbs_keep))+' | cut -f1-3 | sort -V | bedtools merge -i stdin | gzip > '+opts.out
    cmds.append(cmd)
    print '\n'.join(cmds)
    qsub_a_command('qqqq'.join(cmds),opts.out+'_script.sh','qqqq','3G')

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
