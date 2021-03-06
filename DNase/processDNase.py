from optparse import OptionParser
import os
import math
from time import gmtime, strftime
import re
'''
Author:Oana Ursu
Call peaks for Pritchard lab Degner DNase.
'''

def main():
    parser=OptionParser()
    parser.add_option('--dnase_dir',dest='dnase_dir',help='DNase directory. Has /', default='/srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/Degner_mapped_reads/eqtl.uchicago.edu/dsQTL_data/MAPPED_READS/')
    parser.add_option('--peaks_on_subsampled',dest='peaks_on_subsampled',action='store_true')
    parser.add_option('--step_to_perform',dest='step_to_perform',help='liftOver,TrimFromSummit')
    parser.add_option('--out_dir',dest='out_dir',help='',default='/srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/')
    parser.add_option('--genomesize',dest='genomesize',default='2.7e9',help='Genomesize for MACS. Default = 2.7e9 (hs)')
    parser.add_option('--pval_MACS',dest='pval',default='1e-2',help='Pvalue threshold for MACS. Default=1e-2')
    parser.add_option('--subsample_reads',dest='nreads',help='How many million reads')
    parser.add_option('--merge_type',dest='merge_type',help='How to merge peaks: mergeBed, TrimFromSummit, FWHM')
    parser.add_option('--mergedPeaks',dest='mergedPeaks',help='File with merged peaks', default='/srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/peaks/subsampled/mergedPeaks_TrimFromSummit/DNase_TrimFromSummit.mergeBedremoveBlacklist.gz')
    parser.add_option('--total_window',dest='total_window',help='Total window for align2rawsignal. DEFAULT=150',default='150')
    parser.add_option('--seqdir',dest='seqdir',help='Seqdir for align2rawsignal. Default=/srv/gsfs0/projects/kundaje/users/akundaje/projects/encode/data/byDataType/sequence/encodeHg19Male',default='/srv/gsfs0/projects/kundaje/users/akundaje/projects/encode/data/byDataType/sequence/encodeHg19Male')
    parser.add_option('--umap',dest='umap',help='umap directory. DEFAULT=/srv/gsfs0/projects/kundaje/users/akundaje/projects/encode/data/byDataType/umap/hg19_allmappable/globalmap_k1tok1000',default='/srv/gsfs0/projects/kundaje/users/akundaje/projects/encode/data/byDataType/umap/hg19_allmappable/globalmap_k1tok1000')
    parser.add_option('--bigwigdir',dest='bigwigdir',default='/srv/gsfs0/projects/snyder/oursu/histoneQTL/DNase/data/hg19/subsampled/bigwig/')
    opts,args=parser.parse_args()
    os.system('mkdir '+opts.out_dir)

    replicate_groups={}
    for dnase in os.listdir(opts.dnase_dir):
        #All filename things here
        cur_file=opts.dnase_dir+dnase
        out_dir=opts.out_dir+'liftOver/'
        os.system('mkdir '+out_dir)
        cur_bed=out_dir+os.path.basename(cur_file)+'hg18.bed'
        out=out_dir+os.path.basename(cur_file)+'_hg19.tagAlign.gz'
        start_from=out
        shiftReads_out=opts.out_dir+'shiftedReads/'
        os.system('mkdir '+shiftReads_out)
        subsampled_out=opts.out_dir+'subsampled/'
        os.system('mkdir '+subsampled_out)
        #shifted_reads=shiftReads_out+os.path.basename(start_from)+'_75bp_leftOnPosStrand_RightOnNegStrand.tagAlign.bed.gz'
        peaks_dir=opts.out_dir+'peaks/'
        os.system('mkdir '+peaks_dir)
        

        #pooling replicates here
        person=dnase.split('_')[0]
        if person not in replicate_groups.keys():
            replicate_groups[person]=list()
        shifted_reads=shiftReads_out+re.sub('.tagAlign.gz','',os.path.basename(start_from))+'_Rep_'+str(len(replicate_groups[person])+1)+'_'+'.75bp_shifted.tagAlign.gz'
        replicate_groups[person].append(shifted_reads)

        #LIFT OVER
        if opts.step_to_perform=='liftOver':
            cmd_to_bed="zcat "+cur_file+" | awk '{a=$2+20}{print "+'$1"\\t"$2"\\t"a"\\t"'+'"'+os.path.basename(cur_file)+'_read'+'"'+'NR"\\t"1"\\t"$3'+"}' > "+cur_bed
            liftOverData='/home/oursu/devtools/liftOverData/hg18ToHg19.over.chain'
            cmd_liftOver='/home/oursu/devtools/liftOver '+cur_bed+' '+liftOverData+' '+out+'.bed'+' '+out+'_unmapped'
            cmd_zip='cat '+out+'.bed'+' | gzip > '+out
            cmd_rm='rm '+out+'.bed'
            cmds=[cmd_to_bed,cmd_liftOver,cmd_zip,cmd_rm]
            qsub_a_command('qqqq'.join(cmds),out+'_script.sh','qqqq','3G')
        #SHIFT READS
        if opts.step_to_perform=='shiftReads':
            os.system('mkdir '+shiftReads_out)
            shift_dnase_reads(start_from,shifted_reads)
    
    #=========
    #SUBSAMPLE
    #=========
    print replicate_groups
    if opts.step_to_perform=='subsample':
        fileWithReplicates_name=subsampled_out+'DNASE_fileWithReplicates'
        filewithReplicates=open(fileWithReplicates_name,'w')
        for replicate_group in replicate_groups.keys():
            replicates=list(replicate_groups[replicate_group])
            print replicates
            replicates_basename=set()
            for replicate_idx in range(len(replicates)):
                replicate=replicates[replicate_idx]
                replicates_basename.add(os.path.basename(replicate))
            filewithReplicates.write(';'.join(replicates_basename)+'\t'+replicate_group+'\n')
        filewithReplicates.close()
        print fileWithReplicates_name
        #Now, call subsample for all samples
        idir=shiftReads_out[:(len(shiftReads_out)-1)]
        odir=subsampled_out[:(len(subsampled_out)-1)]
        genome='/srv/gs1/projects/kundaje/oursu/Alignment/data/ENCODE_genomes/male/ref.fa.fai'
        nreads=int(float(opts.nreads)*1000000)
        cmd='/srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/subsample_tagAlign_replicates.sh '+idir+' '+odir+' '+fileWithReplicates_name+' '+str(nreads)+' '+genome
        print cmd
        os.system(cmd)
        

    #CALL PEAKS (replicates are pooled)
    # - make pairFile                                                                                                                                              
    # - make fragLenFile             
    # - call MACS2 script from Anshul
    if opts.step_to_perform=='MACS':
        if opts.peaks_on_subsampled:
            print 'Calling peaks on subsampled'
            #Make a subsampled peaks dir
            subsampled_peaks_dir=peaks_dir+'subsampled/'
            os.system('mkdir '+subsampled_peaks_dir)
            pairFileName=subsampled_peaks_dir+'DNASE_pairFile'
            pairFile=open(pairFileName,'w')
            fragLenFileName=subsampled_peaks_dir+'DNASE_fragLenFile'
            fragLenFile=open(fragLenFileName,'w')
            nreads=str(int(float(opts.nreads)*1000000))
            for replicate_group in replicate_groups.keys():
                cur_file=replicate_group+'.subsampleTo'+nreads+'.tagAlign.gz'
                pairFile.write(cur_file+'\n')
                fragLenFile.write(cur_file+'\t'+'150'+'\n')
            pairFile.close()
            fragLenFile.close()
            idir=subsampled_out#[:(len(subsampled_out)-1)]
            odir=subsampled_peaks_dir
            genome='hs'
            chrSizes='/srv/gs1/projects/kundaje/oursu/Alignment/data/ENCODE_genomes/male/ref.fa.fai'
            memory_num='20'
            cmd='/srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/macs2.signal.lsf.submitscript_2014_03_22.sh '+' '+idir+' '+pairFileName+' '+odir+' '+genome+' '+chrSizes+' '+memory_num+' '+fragLenFileName
            print cmd
            os.system(cmd)

        else:
            print 'Calling peaks for everyone'
            pairFile=open(peaks_dir+'DNASE_pairFile','w')
            fragLenFile=open(peaks_dir+'DNASE_fragLenFile','w')
            for replicate_group in replicate_groups.keys():
                replicates=list(replicate_groups[replicate_group])
                current_pairFile_replicates=set()
                for replicate_idx in range(len(replicates)):
                    replicate=replicates[replicate_idx]
                    print os.path.basename(replicate)
                    current_pairFile_replicates.add(os.path.basename(replicate))
                    fragLenFile.write(os.path.basename(replicate)+'\t'+'150'+'\n')
                pairFile.write(';'.join(current_pairFile_replicates)+'\n')
            pairFile.close()
            fragLenFile.close()
            idir=shiftReads_out
            pairFileName=peaks_dir+'DNASE_pairFile'
            fragLenFileName=peaks_dir+'DNASE_fragLenFile'
            odir=peaks_dir
            genome='hs'
            chrSizes='/srv/gs1/projects/kundaje/oursu/Alignment/data/ENCODE_genomes/male/ref.fa.fai'
            memory_num='20'
            cmd='/srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/macs2.signal.lsf.submitscript_2014_03_22.sh '+' '+idir+' '+pairFileName+' '+odir+' '+genome+' '+chrSizes+' '+memory_num+' '+fragLenFileName
            print cmd
            os.system(cmd)

    if opts.step_to_perform=='TrimFromSummit':
        subsampled_peaks_dir=peaks_dir+'subsampled/'
        for replicate_group in replicate_groups.keys():
            nreads=str(int(float(opts.nreads)*1000000))
            cmd0='module load bedtools/2.18.0'
            chrSizes='/srv/gs1/projects/kundaje/oursu/Alignment/data/ENCODE_genomes/male/ref.fa.fai'
            cur_file=subsampled_peaks_dir+replicate_group+'.subsampleTo'+nreads+'.tagAlign'+'_VS_NoControl_peaks.narrowPeak'
            trimCmd="cat "+cur_file+' | awk \'{trimstart=$2+$10-50} {trimend=$2+$10+50} {print $1"\\t"trimstart"\\t"trimend"\\t"$8}\' | bedtools slop -b 0 -i stdin -g '+chrSizes+' > '+cur_file+'.100bpCenterSummit.bed'
            qsub_a_command('qqqq'.join([cmd0,trimCmd]),cur_file+'trimCmd_script.sh',split_string='qqqq',memory_number='3G')

    if opts.step_to_perform=='extractSignalPerIndividual_trimmedFromSummit':
        subsampled_peaks_dir=peaks_dir+'subsampled/'
        for replicate_group in replicate_groups.keys():
            nreads=str(int(opts.nreads)*1000000)
            replicate_col1=replicate_group+'.subsampleTo'+nreads+'.tagAlign.gz'
            replicate_col2=re.sub('tagAlign.gz','wiggler',replicate_col1)
            align2rawsignal_mat=subsampled_out+'/align2rawsignal/'+replicate_col2+'.norm5.rawsignal.mat'
            cur_file=subsampled_peaks_dir+replicate_group+'.subsampleTo'+nreads+'.tagAlign'+'_VS_NoControl_peaks.narrowPeak'+'.100bpCenterSummit.bed'
            cmds=list()
            output_extract_dir=subsampled_out+'/align2rawsignal/extractSignal_peaksOfEachIndividual/'
            os.system('mkdir '+output_extract_dir)
            output_extract=output_extract_dir+os.path.basename(align2rawsignal_mat)+'_VS_'+replicate_group+'_peaks100bpCenterSummit.MEAN.cagt'
            chrSizes='/srv/gs1/projects/kundaje/oursu/Alignment/data/ENCODE_genomes/male/ref.fa.fai'
            cmds.append('module load bedtools/2.18.0')
            #cmds.append('cat '+cur_file+' | cut -f1-3 | bedtools slop -b 0 -i stdin -g '+chrSizes+' > '+cur_file+'.bed')
            cmds.append('source  /srv/gsfs0/projects/kundaje/commonRepository/src/lab_bashrc')
            cmds.append('activateMCR')
            cmds.append('extractSignal -i='+cur_file+' -t='+align2rawsignal_mat+' -mf=mean -of=cagt -if=bed -o='+output_extract)
            cmds.append('cat '+output_extract+' | sort -rn -k4 | head -n1000000 > '+output_extract+'.top1M.bed')
            print '\n'.join(cmds)
            qsub_a_command('qqqq'.join(cmds),output_extract_dir+replicate_group+'_script.sh','qqqq',memory_number='40G')
            

    #extractSignal
    if opts.step_to_perform=='extractSignal':
        cur_merged_peaks=opts.mergedPeaks
        first=True
        for replicate_group in replicate_groups.keys():
            nreads=str(int(opts.nreads)*1000000)
            replicate_col1=replicate_group+'.subsampleTo'+nreads+'.tagAlign.gz'
            replicate_col2=re.sub('tagAlign.gz','wiggler',replicate_col1)
            align2rawsignal_mat=subsampled_out+'/align2rawsignal/'+replicate_col2+'.norm5.rawsignal.mat'
            cmds=list()
            output_extract_dir=subsampled_out+'/align2rawsignal/extractSignal/'
            os.system('mkdir '+output_extract_dir)
            output_extract=output_extract_dir+os.path.basename(align2rawsignal_mat)+'_VS_'+'mergePeaks_'+re.sub('.gz','',os.path.basename(cur_merged_peaks))+'.MEAN.cagt'
            cmds.append('source  /srv/gsfs0/projects/kundaje/commonRepository/src/lab_bashrc')
            cmds.append('activateMCR')
            if first:
                cmds.append('zcat '+cur_merged_peaks+' > '+cur_merged_peaks+'_notgz')
                first=False
            cmds.append('extractSignal -i='+cur_merged_peaks+'_notgz'+' -t='+align2rawsignal_mat+' -mf=mean -of=cagt -if=bed -o='+output_extract)                        
            #cmds.append('rm '+cur_merged_peaks+'_notgz')
            print '\n'.join(cmds)
            qsub_a_command('qqqq'.join(cmds),output_extract_dir+replicate_group+'_script.sh','qqqq',memory_number='40G')

    if opts.step_to_perform=='align2rawsignal':
        seqdir=opts.seqdir
        umap=opts.umap
        indir=subsampled_out
        for replicate_group in replicate_groups.keys():
            print replicate_group
            outdir=indir+'/align2rawsignal'
            os.system('mkdir '+outdir)
            nreads=str(int(float(opts.nreads)*1000000))
            replicate_col1=replicate_group+'.subsampleTo'+nreads+'.tagAlign.gz'
            replicate_col2=re.sub('tagAlign.gz','wiggler',replicate_col1)
            cur_pairfile=outdir+'/'+replicate_col2+'.pairFile'
            os.system('echo "'+replicate_col1+'\t'+replicate_col2+'" > '+cur_pairfile)
            cur_fragLen=outdir+'/'+replicate_col2+'.fragLen'
            total_window='150'
            cur_fragLen_value='150'
            os.system('echo "'+replicate_col1+'\t'+cur_fragLen_value+'\t'+total_window+'" > '+cur_fragLen)
            export_tmp_cmd='export TMP='+outdir+replicate_col2+'TMP'
            align2rawsignal_cmd='/srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/wiggler/batchWiggler.sge.sh '+cur_pairfile+' '+cur_fragLen+' '+indir+' '+outdir+' TRUE 20 mat epanechnikov 5 0 '+seqdir+' '+umap
            os.system(export_tmp_cmd)
            print align2rawsignal_cmd
            os.system(align2rawsignal_cmd)

    if opts.step_to_perform=='bigwig':
        seqdir=opts.seqdir
        umap=opts.umap
        indir=subsampled_out
        for replicate_group in replicate_groups.keys():
            print replicate_group
            outdir=indir+'/align2rawsignal'
            os.system('mkdir '+outdir)
            nreads=str(int(float(opts.nreads)*1000000))
            replicate_col1=replicate_group+'.subsampleTo'+nreads+'.tagAlign.gz'
            replicate_col2=re.sub('tagAlign.gz','wiggler',replicate_col1)
            cur_pairfile=outdir+'/'+replicate_col2+'.pairFile'
            os.system('echo "'+replicate_col1+'\t'+replicate_col2+'" > '+cur_pairfile)
            cur_fragLen=outdir+'/'+replicate_col2+'.fragLen'
            total_window='150'
            cur_fragLen_value='150'
            os.system('echo "'+replicate_col1+'\t'+cur_fragLen_value+'\t'+total_window+'" > '+cur_fragLen)
            export_tmp_cmd='export TMP=/srv/gsfs0/scratch/oursu/'
            chrSizes='/srv/gs1/projects/kundaje/oursu/Alignment/data/ENCODE_genomes/male/ref.fa.fai'
            align2rawsignal_cmd='/srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/wiggler/batchWiggler.sge.sh '+cur_pairfile+' '+cur_fragLen+' '+indir+' '+opts.bigwigdir+' TRUE 20 bg epanechnikov 5 0 '+seqdir+' '+umap+' '+chrSizes
            os.system(export_tmp_cmd)
            print align2rawsignal_cmd
            os.system(align2rawsignal_cmd)



    #Merging options for DNase - merge directly just the peaks that have a good pvalue in at least 1 individual
    if opts.step_to_perform=='mergePeaks':
        merge_peaks=set()
        merge_peaks_TrimFromSummit=set()
        blacklist_file='/srv/gs1/projects/kundaje/oursu/Alignment/data/blacklist/wgEncodeDacMapabilityConsensusExcludable.bed.gz'
        #TODO: provide support for non-subsampled
        if opts.peaks_on_subsampled:
            print 'Merging peaks on subsampled'
            subsampled_peaks_dir=peaks_dir+'subsampled/'
            nreads=str(int(float(opts.nreads)*1000000))
            for replicate_group in replicate_groups.keys():
                cur_file=replicate_group+'.subsampleTo'+nreads+'.tagAlign'+'_VS_NoControl_peaks.narrowPeak'
                merge_peaks.add(subsampled_peaks_dir+cur_file)
                merge_peaks_TrimFromSummit.add(subsampled_peaks_dir+cur_file+'.100bpCenterSummit.bed')
            if opts.merge_type=='TrimFromSummit':
                merged_dir=subsampled_peaks_dir+'mergedPeaks_TrimFromSummit/'
                os.system('mkdir '+merged_dir)
                merge_name=merged_dir+'DNase_TrimFromSummit.mergeBed.gz'
                cmd0='module load bedtools/2.18.0'
                cmd1='cat '+' '.join(list(merge_peaks_TrimFromSummit))+' | sort -V | bedtools merge -i stdin -nms | gzip > '+merge_name
                blacklist_site=blacklist_file+'notgz'
                ungzblacklist='zcat '+blacklist_file+' > '+blacklist_site
                merged_noBL=re.sub('.gz','',merge_name)+'removeBlacklist.gz'
                blacklist_cmd='zcat '+merge_name+' | subtractBed -A -a stdin -b '+blacklist_site+' | gzip > '+merged_noBL
                cmd_rm='rm '+blacklist_site
                #And now, remove peaks with bad quality across all individuals                                                                                           
                maxMinPvalue_minusLog10='5'                                                                                                                                
                cmd_pval='/srv/gs1/software/R/R-3.1.0/bin/Rscript /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/chromoVariationCode/keepPeaks_by_maxMinPvalue.R '+merged_noBL+' '+maxMinPvalue_minusLog10+' '+re.sub('.gz','',merged_noBL)+'_Log10PvalueThreshold_'+maxMinPvalue_minusLog10+'.gz'
                qsub_a_command('qqqq'.join([cmd0,cmd1,ungzblacklist,blacklist_cmd,cmd_rm,cmd_pval]),merge_name+'_script.sh',split_string='qqqq',memory_number='3G')

            if opts.merge_type=='mergeBed':
                merged_dir=subsampled_peaks_dir+'mergedPeaks/'
                os.system('mkdir '+merged_dir)
                merge_name=merged_dir+'DNase_mergeBed.gz'
                cmd0='module load bedtools/2.18.0'
                cmd1='cat '+' '.join(list(merge_peaks))+' | cut -f1-3,8 | sort -V | bedtools merge -i stdin -nms | gzip > '+merge_name
                blacklist_site=blacklist_file+'notgz'
                ungzblacklist='zcat '+blacklist_file+' > '+blacklist_site
                blacklist_cmd='zcat '+merge_name+' | subtractBed -A -a stdin -b '+blacklist_site+' | gzip > '+merge_name+'removeBlacklist.gz'
                cmd_rm='rm '+blacklist_site
                qsub_a_command('qqqq'.join([cmd0,cmd1,ungzblacklist,blacklist_cmd,cmd_rm]),merge_name+'_script.sh',split_string='qqqq',memory_number='3G')

            if opts.merge_type=='mergeBed_topInEachIndividual_byCol':
                merged_dir=subsampled_peaks_dir+'mergedPeaks_TrimFromSummit_top1MSignalsPerIndividual/'
                os.system('mkdir '+merged_dir)
                merge_name=merged_dir+'DNase_TrimFromSummit_top1MSignalsPerIndividual.mergeBed.gz'
                cmd0='module load bedtools/2.18.0'
                #make a file with top things from each person
                cmd1='cat '+' '.join(list(merge_peaks_TrimFromSummit))+' | sort -V | bedtools merge -i stdin -nms | gzip > '+merge_name
                blacklist_site=blacklist_file+'notgz'
                ungzblacklist='zcat '+blacklist_file+' > '+blacklist_site
                merged_noBL=re.sub('.gz','',merge_name)+'removeBlacklist.gz'
                blacklist_cmd='zcat '+merge_name+' | subtractBed -A -a stdin -b '+blacklist_site+' | gzip > '+merged_noBL
                cmd_rm='rm '+blacklist_site



def shift_dnase_reads(start_from,shifted_reads):
    genome_size='/srv/gs1/projects/kundaje/oursu/Alignment/data/ENCODE_genomes/male/ref.fa.fai'
    #Make 2 files, one with + strand one with -
    print 'Shifting reads for '+start_from
    plus_strand=start_from+'plusStrand'
    minus_strand=start_from+'minusStrand'
    minus_strand_slop=minus_strand+'_slopBed'
    plus_strand_slop=plus_strand+'_slopBed'
    cmd_plus='zcat '+start_from+' | grep + > '+plus_strand
    cmd_minus='zcat '+start_from+' | grep - > '+minus_strand
    #On +, subtract 75
    cmd_shift_plus='/srv/gs1/software/bedtools/bedtools-2.17.0/bin/slopBed -l 75 -r -75 -i '+plus_strand+' -g '+genome_size+' > '+plus_strand_slop
    #On -, add 75
    cmd_shift_minus='/srv/gs1/software/bedtools/bedtools-2.17.0/bin/slopBed -l -75 -r 75 -i '+minus_strand+' -g '+genome_size+' > '+minus_strand_slop
    chr_grep="-w 'chr1\|chr2\|chr3\|chr4\|chr5\|chr6\|chr7\|chr8\|chr9\|chr10\|chr11\|chr12\|chr13\|chr14\|chr15\|chr16\|chr17\|chr18\|chr19\|chr20\|chr21\|chr22\|chrX\|chrY'"
    combine_cmd='cat '+minus_strand_slop+' '+plus_strand_slop+' | grep '+chr_grep+' | /home/oursu/devtools/bedClip stdin '+genome_size+' '+shifted_reads+'clipped'
    final_cmd='cat '+shifted_reads+'clipped'+' | gzip > '+shifted_reads
    rm_plus='rm '+plus_strand
    rm_minus='rm '+minus_strand
    rm_plus_slop='rm '+plus_strand_slop
    rm_minus_slop='rm '+minus_strand_slop
    rm_clipped='rm '+shifted_reads+'clipped'
    cmds=[cmd_plus,cmd_minus,cmd_shift_plus,cmd_shift_minus,combine_cmd,final_cmd,rm_plus,rm_minus,rm_plus_slop,rm_minus_slop,rm_clipped]
    print '\n'.join(cmds)
    qsub_a_command('qqqq'.join(cmds),shifted_reads+'_script.sh','qqqq','3G')

def qsub_a_command(cmd,shell_script_name,split_string=',',memory_number='20G'):
    f=open(shell_script_name,'w')
    cmds=cmd.split(split_string)
    for i in range(len(cmds)):
        f.write(cmds[i]+'\n') #just write the command
    f.close()
    #make runnable                                                                                                                                                      
    os.system('chmod 711 '+shell_script_name)
    #Qsub the script                                                                                                                                                    
    os.system("qsub -l mem_free="+memory_number+" -l h_vmem="+memory_number+" -l h_rt=20:00:00 -o "+shell_script_name+'.o'+' -e '+shell_script_name+'.e'+' '+shell_script_name)
        

main()
