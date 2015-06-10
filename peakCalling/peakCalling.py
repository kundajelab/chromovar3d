from optparse import OptionParser
import os
import math
from time import gmtime, strftime
import re
import sys
'''
Author:Oana Ursu
Call peaks.
'''

def main():
    parser=OptionParser()
    
    parser.add_option('--code_path',dest='code_path',help='Path of the code, to make it easy to transfer code and have it still work. DEFAULT: /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/', default='/srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/')
    parser.add_option('--metadata',dest='metadata',help='Metadata file. One line per condition. Should be tab or space-delimited: 1. Individual, 2. sample name (unique),3. fastq1, 4. fastq2, 5. genome_path (for instance for <path>/NA19099 the genome_path=path), 6. gender,7. vcf file for personal genome,alignment directory. If any of these entries is missing, e.g. fastq2 is missing, say NA. Header should start with #. NOTE: input file must be in the same directory as the samples analyzed!!!')
    parser.add_option('--columns_for_chip',dest='cols_chip',help='Columns of interest from the metadata file for Chip peak calling.<sample_name>,<bam>,<qc_dir>,<peak_dir>,<replicate group>,<isinput>,<inputName>,<mergeGroup>,<mergeDir>. Columns start at 0. Default: 1,9,10,11,13,14,15,16,17.',default='1,9,10,12,13,14,15,16,17')
    parser.add_option('--step_to_perform',dest='step_to_perform',help='Step to perform. options: SPP, mergeReplicates, subsample, MACS, align2rawsignal, align2rawsignal_bigwig,mergePeaks, extractsignal.')
    parser.add_option('--sample_names_to_do',dest='todo',help='Sample names for subset of things to run',default='')
    parser.add_option('--genomesize',dest='genomesize',default='hs',help='Genomesize for MACS. Default = hs')
    parser.add_option('--chrSizes',dest='chrSizes',default='/srv/gs1/projects/kundaje/oursu/Alignment/data/ENCODE_genomes/male/ref.fa.fai',help='ChrSizes default=/srv/gs1/projects/kundaje/oursu/Alignment/data/ENCODE_genomes/male/ref.fa.fai')
    parser.add_option('--pval_MACS',dest='pval',default='1e-2',help='Pvalue threshold for MACS. Default=1e-2')
    #parser.add_option('--merge_dir',dest='merge_dir',default='',help='Directory to put merged peaks. default=directory for peaks of first entry in metadata')
    parser.add_option('--subsample_reads',dest='subsample_reads',help='Only matters if you subsample reads. Then you specify millions of reads desired. Else, it will be treated as empty')
    parser.add_option('--peaks_on_subsampled',dest='peaks_on_subsampled',action='store_true',help='Call peaks on the subsampled version. Will give an error if --subsample_reads is not given')
    parser.add_option('--blacklist',dest='blacklist',help='Provide a blacklist if you want your peaks to be filtered for the blacklist. ',default='/srv/gs1/projects/kundaje/oursu/Alignment/data/blacklist/wgEncodeDacMapabilityConsensusExcludable.bed.gz')
    parser.add_option('--total_window',dest='total_window',help='Total window for align2rawsignal. DEFAULT=150',default='150')
    parser.add_option('--seqdir',dest='seqdir',help='Seqdir for align2rawsignal. Default=/srv/gsfs0/projects/kundaje/users/akundaje/projects/encode/data/byDataType/sequence/encodeHg19Male',default='/srv/gsfs0/projects/kundaje/users/akundaje/projects/encode/data/byDataType/sequence/encodeHg19Male')
    parser.add_option('--umap',dest='umap',help='umap directory. DEFAULT=/srv/gsfs0/projects/kundaje/users/akundaje/projects/encode/data/byDataType/umap/hg19_allmappable/globalmap_k1tok1000',default='/srv/gsfs0/projects/kundaje/users/akundaje/projects/encode/data/byDataType/umap/hg19_allmappable/globalmap_k1tok1000')
    parser.add_option('--align2rawsignal_bigwigdir',dest='align2rawsignal_bigwigdir',default='/srv/gsfs0/projects/snyder/oursu/histoneQTL/ChIPseq_alignment/results/peakCalls/align2rawsignal_bigwig')
    opts,args=parser.parse_args()
    

    [sample_name_col_str,bam_col_str,qc_dir_col_str,peak_dir_col_str,replicate_col_str,input_col_str,input_name_str,merge_group_str,merge_dir_str]=opts.cols_chip.split(',')
    bam_col=int(bam_col_str)
    qc_dir_col=int(qc_dir_col_str)
    peak_dir_col=int(peak_dir_col_str)
    sample_name_col=int(sample_name_col_str)
    replicate_col=int(replicate_col_str)
    input_col=int(input_col_str)
    input_name_col=int(input_name_str)
    merge_group_col=int(merge_group_str)
    merge_dir_col=int(merge_dir_str)

    sample_di={}
    for line in open(opts.metadata).readlines():
        if line[0]=='#':
            continue
        items=line.strip().split()
        sample_name=items[sample_name_col]
        bam=items[bam_col]
        qc_dir=items[qc_dir_col]
        peak_dir=items[peak_dir_col]
        sample_di[sample_name]={}
        sample_di[sample_name]['bam']=bam
        sample_di[sample_name]['QC_dir']=qc_dir
        sample_di[sample_name]['peak_dir']=peak_dir
        sample_di[sample_name]['replicate_group']=items[replicate_col]
        sample_di[sample_name]['input']=items[input_name_col]
        sample_di[sample_name]['merge_group']=items[merge_group_col]
        sample_di[sample_name]['merge_dir']=items[merge_dir_col]
    if opts.todo=='':
        of_interest=sample_di.keys()
    if opts.todo!='':
        of_interest_items=opts.todo.split(',')
        of_interest=set()
        for sample_name in sample_di.keys():
            for of_interest_item in of_interest_items:
                if of_interest_item in sample_name:
                    of_interest.add(sample_name)
    #=====
    # SPP
    #=====
    samtools_setup_cmd='export PATH=/srv/gs1/software/samtools/samtools-0.1.19/bin/:$PATH'
    if opts.step_to_perform=='SPP':
        for sample_name in sample_di.keys():
            if sample_name not in of_interest:
                continue
            spp_cmd='/home/oursu/devtools/R-3.0.2/bin/Rscript /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/phantompeakqualtools/run_spp_nodups.R -c='+sample_di[sample_name]['bam']+' '+'-odir='+sample_di[sample_name]['QC_dir']+' '+'-savp -out='+sample_di[sample_name]['QC_dir']+os.path.basename(sample_di[sample_name]['bam'])+'SPP_table.txt'
            cmds=samtools_setup_cmd+'qqqq'+spp_cmd
            print cmds
            qsub_a_command(cmds,sample_di[sample_name]['QC_dir']+os.path.basename(sample_di[sample_name]['bam'])+'SPP_script.sh','qqqq','3G')
    if opts.step_to_perform=='SPP':
        sys.exit()
    #put replicates into replicate groups and pool all inputs together                                                                                                   
    replicate_groups={}
    for sample_name in sample_di.keys():
        if sample_name not in of_interest:
            continue
        replicate_group=sample_di[sample_name]['replicate_group']
        if replicate_group not in replicate_groups.keys():
            replicate_groups[replicate_group]={}
            #replicate_groups[replicate_group]['fragLens']=set()
            replicate_groups[replicate_group]['bam']=set()
            replicate_groups[replicate_group]['peak_dir']=sample_di[sample_name]['peak_dir']
            replicate_groups[replicate_group]['input_name']=set()
            replicate_groups[replicate_group]['merge_group']=sample_di[sample_name]['merge_group']
            replicate_groups[replicate_group]['merge_dir']=sample_di[sample_name]['merge_dir']
        replicate_groups[replicate_group]['bam'].add(sample_di[sample_name]['bam'])
        replicate_groups[replicate_group]['input_name'].add(sample_di[sample_name]['input'])

    #=====================
    #MERGE REPLICATES HERE
    #=====================
    #Put all replicates into 1 file. Then run SPP on that file. Take the output from SPP and create a new fragLen file for everyone.
    spp_tables=set()
    #Setup directory for all merged replicates
    dir_replicate0=os.path.dirname(list(replicate_groups[replicate_groups.keys()[0]]['bam'])[0])+'/merged_replicates/'
    os.system('mkdir '+dir_replicate0)
    os.system('mkdir '+dir_replicate0+'QC')
    if opts.step_to_perform=='mergeReplicates':
        spp_tables=set()
        for replicate_group in replicate_groups.keys():
            merged_replicate_name=dir_replicate0+'Rep_0_'+replicate_group+'.mergedReplicates.tagAlign.gz'
            curr_spp_table=dir_replicate0+'QC/'+os.path.basename(merged_replicate_name)+'SPP_table.txt'
            spp_tables.add(curr_spp_table)
            things_to_merge=re.sub('.bam','.tagAlign.gz',' '.join(replicate_groups[replicate_group]['bam']))
            cmd_mergeReplicates='zcat '+things_to_merge+' | gzip > '+merged_replicate_name
            #cmd_spp_mergeReplicates='/home/oursu/devtools/R-3.0.2/bin/Rscript /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/phantompeakqualtools/run_spp_nodups.R -c='+merged_replicate_name+' '+'-odir='+dir_replicate0+'QC'+' '+'-savp -out='+curr_spp_table
            cmd_spp_mergeReplicates='' #For now, don't spp merged replicates, we will spp the subsampled things anyway
            print cmd_mergeReplicates
            print cmd_spp_mergeReplicates
            qsub_a_command(cmd_mergeReplicates+'qqqq'+cmd_spp_mergeReplicates,merged_replicate_name+'_script.sh','qqqq','3G')

    #=========
    #SUBSAMPLE
    #=========
    subsample_dir=dir_replicate0+'subsampled/'
    os.system('mkdir '+subsample_dir)
    os.system('mkdir '+subsample_dir+'QC')
    if opts.step_to_perform=='subsample':
        print 'subsampling'
        nreads=opts.subsample_reads
        if nreads=='':
            sys.exit('You want to subsample but you have not specified how many reads you desire. Please specify the subsample_reads parameter')
        for replicate_group in replicate_groups.keys():
            merged_replicate_name=dir_replicate0+'Rep_0_'+replicate_group+'.mergedReplicates.tagAlign.gz'
            #1st, from merged replicate figure out the percent reads we need
            cmds=[]
            n_reads_given=str(int(nreads)*1000000)
            cmds.append('module load samtools/0.1.19')
            cmds.append('module add bedtools/2.18.0')
            cmds.append('seed_value=7')
            get_total_reads='OREADS=$(zcat '+merged_replicate_name+' | wc -l)'
            cmds.append('NREADS='+n_reads_given)
            cmds.append(get_total_reads)
            get_proportion='PRO=$(echo "scale=6;$NREADS/$OREADS" | bc)'
            cmds.append(get_proportion)
            #For bam in each replicate group, subsample and convert to tagAlign directly. Then merge them all into 1.
            result_files=[]
            reads_per_bam=int(math.ceil(int(n_reads_given)/len(replicate_groups[replicate_group]['bam'])))
            for bam in replicate_groups[replicate_group]['bam']:
                result_file=subsample_dir+os.path.basename(re.sub('.bam','',bam)+'.subsampleTo'+n_reads_given+'.tagAlign.gz')
                cmds.append('S_USED=""')
                cmds.append('if (( $(echo "$PRO < 1.0" | bc -l) )); then S_USED=" -s ${seed_value}${PRO}"; fi')
                cmds.append("samtools view${S_USED} -b "+bam+" | bamToBed -i stdin | awk 'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{$5=\"1000\" ; print $0}' | gzip > "+result_file)
                result_files.append(result_file)
            #Merge the subsampled tagAligns
            output_name=subsample_dir+re.sub('.tagAlign.gz','',os.path.basename(merged_replicate_name))+'.subsampleTo'+n_reads_given+'.tagAlign.gz'
            cmds.append('zcat '+' '.join(result_files)+' | gzip > '+output_name)
            curr_spp_table=subsample_dir+'QC/'+os.path.basename(output_name)+'SPP_table.txt'
            cmds.append('/home/oursu/devtools/R-3.0.2/bin/Rscript /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/phantompeakqualtools/run_spp_nodups.R -c='+output_name+' '+'-odir='+subsample_dir+'QC'+' '+'-savp -out='+curr_spp_table)
            print '\n'.join(cmds)
            qsub_a_command('qqqq'.join(cmds),output_name+'_script.sh','qqqq','3G')
            print '---'
            
    #========================================================================
    # MACS. 2 options. Either run the subsampled reads or the original files.
    #========================================================================
    #Run Anshul's code
    load_macs_cmd='module load MACS2/2.0.10'
    if opts.step_to_perform=='MACS':
        if opts.peaks_on_subsampled:
            n_reads_given=str(int(opts.subsample_reads)*1000000)
            print 'Calling peaks on subsampled'
            idir=subsample_dir
            genome=opts.genomesize
            chrsizes=opts.chrSizes
            memo='8'
            for replicate_group in replicate_groups.keys():
                odir=replicate_groups[replicate_group]['peak_dir']
                os.system('mkdir '+odir)
                replicate_col1=re.sub('.tagAlign.gz','','Rep_0_'+replicate_group+'.mergedReplicates.tagAlign.gz')+'.subsampleTo'+n_reads_given+'.tagAlign.gz'
                cur_input=';'.join(replicate_groups[replicate_group]['input_name'])
                print "Removing directory information for the input name, assuming that the input file is in the directory idir. TODO: perhaps change this for people who don't want this functionality"
                cur_input=re.sub(os.path.dirname(list(replicate_groups[replicate_group]['input_name'])[0])+'/','',cur_input)
                cur_pairfile=odir+replicate_col1+'.pairFile'
                os.system('echo "'+replicate_col1+'\t'+cur_input+'" > '+cur_pairfile)
                cur_fragLen=odir+replicate_col1+'.fragLen'
                cur_fragLen_value=open(idir+'/QC/'+replicate_col1+'SPP_table.txt','r').readlines()[0].strip().split()[2].split(',')[0]
                os.system('echo "'+replicate_col1+'\t'+cur_fragLen_value+'" > '+cur_fragLen)
                macs_cmd='/srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/macs2.signal.lsf.submitscript_2014_03_22.sh '+idir+' '+' '+cur_pairfile+' '+odir[0:(len(odir)-1)]+' '+genome+' '+chrsizes+' '+memo+' '+cur_fragLen
                print macs_cmd
                os.system(macs_cmd)
        else:
            print 'Calling peaks on full data. Not yet implemented'


    if opts.step_to_perform=='extractSignal':
        for replicate_group in replicate_groups.keys():
            maxMinPvalue_minusLog10='5'
            if opts.peaks_on_subsampled:
                blacklist_add=''
                if opts.blacklist!='':
                    blacklist_add='removeBlacklist.gz'
                #signal extracted for everything but the datamatrix will be just with blacklist removed peaks
                #Determine the merge file
                #cur_merged_peaks=replicate_groups[replicate_group]['merge_dir']+'mergePeaks_'+replicate_groups[replicate_group]['merge_group']+'.gz'+blacklist_add
                merge_name=replicate_groups[replicate_group]['merge_dir']+'mergePeaks_'+replicate_groups[replicate_group]['merge_group']+'.gz'
                cur_merged_peaks=merge_name+'removeBlacklist.gz'+'_Log10PvalueThreshold_'+maxMinPvalue_minusLog10+'.gz'
                n_reads_given=str(int(opts.subsample_reads)*1000000)
                replicate_col1=re.sub('.tagAlign.gz','','Rep_0_'+replicate_group+'.mergedReplicates.tagAlign.gz')+'.subsampleTo'+n_reads_given+'.tagAlign.gz'
                replicate_col2=re.sub('tagAlign.gz','wiggler',replicate_col1)
                align2rawsignal_mat=replicate_groups[replicate_group]['peak_dir']+'align2rawsignal/'+replicate_col2+'.norm5.rawsignal.mat'
                cmds=list()
                #output_extract_dir=replicate_groups[replicate_group]['peak_dir']+'align2rawsignal/extractSignal'+blacklist_add+'/'
                output_extract_dir=replicate_groups[replicate_group]['merge_dir']+'extractSignal/'
                os.system('mkdir '+output_extract_dir)
                output_extract=output_extract_dir+os.path.basename(align2rawsignal_mat)+'_VS_'+'mergePeaks_'+os.path.basename(cur_merged_peaks)+'.MEAN.cagt'
                cmds.append('source  /srv/gsfs0/projects/kundaje/commonRepository/src/lab_bashrc')
                cmds.append('activateMCR')
                cmds.append('zcat '+cur_merged_peaks+' > '+cur_merged_peaks+'_notgz') 
                cmds.append('extractSignal -i='+cur_merged_peaks+'_notgz'+' -t='+align2rawsignal_mat+' -mf=mean -of=cagt -if=bed -o='+output_extract)
                cmds.append('rm '+cur_merged_peaks+'_notgz')
                cmds.append('cat '+output_extract+' | gzip > '+output_extract+'.gz')
                cmds.append('rm '+output_extract)
                print '\n'.join(cmds)
                qsub_a_command('qqqq'.join(cmds),output_extract_dir+replicate_group+'_script.sh','qqqq',memory_number='20G')
                

    #Done with /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/chromoVariationCode/mergeSignalValues.sh
    '''
    if opts.step_to_perform=='combineSignal':
        blacklist_add=''
        if opts.blacklist!='':
            blacklist_add='removeBlacklist.gz' 
        data_dir=replicate_groups[replicate_group]['peak_dir']+'data_matrix/'
        os.system('mkdir '+data_dir)
        merge_replicates={}
        for replicate_group in replicate_groups.keys():
            if replicate_groups[replicate_group]['merge_group'] not in merge_replicates.keys():
                merge_replicates[replicate_groups[replicate_group]['merge_group']]=set()
            merge_replicates[replicate_groups[replicate_group]['merge_group']].add(replicate_group)
        datastart=replicate_groups[replicate_group]['peak_dir']+'align2rawsignal/extractSignal/Rep_0_'
        for merge_group in merge_replicates.keys():
            cur_out=data_dir+merge_group+blacklist_add.strip('.gz')+'_signal'
            if opts.peaks_on_subsampled:
                print 'merge_group'
                print merge_group
                n_reads_given=str(int(opts.subsample_reads)*1000000)
                dataend='.mergedReplicates.subsampleTo'+n_reads_given+'.wiggler.norm5.rawsignal.mat_VS_mergePeaks_'+merge_group+'.MEAN.cagt'
                cmd='python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/merge_signal_values.py --out '+cur_out+' --mark '+merge_group+' --replicate_groups '+','.join(list(merge_replicates[merge_group]))+' --datastart '+datastart+' --dataend '+dataend+' --peakkeep '+replicate_groups[replicate_group]['peak_dir']+'merged_peaks_'+merge_group+'/'+'mergePeaks_'+merge_group+'.gz'+blacklist_add
                qsub_a_command(cmd,data_dir+merge_group+blacklist_add+'_script.sh','qqqq',memory_number='20G')
                print cmd
    '''            

    #===============
    #ALIGN2RAWSIGNAL
    #===============
    if opts.step_to_perform=='align2rawsignal':
        cmds_initial=initialize_wiggler()
        for cmd in cmds_initial:
            os.system(cmd)
        seqdir=opts.seqdir
        umap=opts.umap
        if opts.peaks_on_subsampled:
            print 'Generating signal on subsampled'
            indir=subsample_dir
        else:
            indir=dir_replicate0
        for replicate_group in replicate_groups.keys():
            try:
                if opts.peaks_on_subsampled:
                    outdir=replicate_groups[replicate_group]['peak_dir']+'align2rawsignal'
                    os.system('mkdir '+outdir)
                    n_reads_given=str(int(opts.subsample_reads)*1000000)
                    replicate_col1=re.sub('.tagAlign.gz','','Rep_0_'+replicate_group+'.mergedReplicates.tagAlign.gz')+'.subsampleTo'+n_reads_given+'.tagAlign.gz'
                    replicate_col2=re.sub('tagAlign.gz','wiggler',replicate_col1)
                    cur_pairfile=outdir+'/'+replicate_col2+'.pairFile'
                    os.system('echo "'+replicate_col1+'\t'+replicate_col2+'" > '+cur_pairfile)
                    cur_fragLen=outdir+'/'+replicate_col2+'.fragLen'
                    total_window=opts.total_window
                    cur_fragLen_value=open(indir+'/QC/'+replicate_col1+'SPP_table.txt','r').readlines()[0].strip().split()[2].split(',')[0]
                    os.system('echo "'+replicate_col1+'\t'+cur_fragLen_value+'\t'+total_window+'" > '+cur_fragLen)
                    export_tmp_cmd='export TMP='+outdir+replicate_col2+'TMP'
                    align2rawsignal_cmd='/srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/wiggler/batchWiggler.sge.sh '+cur_pairfile+' '+cur_fragLen+' '+indir+' '+outdir+' TRUE 20 mat epanechnikov 5 0 '+seqdir+' '+umap+' '+opts.chrSizes
                    os.system(export_tmp_cmd)
                    print align2rawsignal_cmd
                    os.system(align2rawsignal_cmd)
                    
                    '''
                    #print '\n'.join(cmds_initial)+'\n'+export_tmp_cmd+'\n'+align2rawsignal_cmd
                    #qsub_a_command('qqqq'.join(cmds_initial)+'qqqq'+export_tmp_cmd+'qqqq'+align2rawsignal_cmd,outdir+replicate_col2+'_script.sh','qqqq',memory_number='2G')
                    '''
                else:
                    print 'Under construction to call on the not subsampled things'
            except:
                print "not running "+replicate_group

    if opts.step_to_perform=='bigwig':
        cmds_initial=initialize_wiggler()
        for cmd in cmds_initial:
            os.system(cmd)
        seqdir=opts.seqdir
        umap=opts.umap
        if opts.peaks_on_subsampled:
            print 'Generating signal on subsampled'
            indir=subsample_dir
        else:
            indir=dir_replicate0
        for replicate_group in replicate_groups.keys():
            print 'new----------------------------------------'
            print 'subsample dir'
            print subsample_dir
            print 'repl 0 '
            print dir_replicate0
            try:
                if opts.peaks_on_subsampled:
                    outdir=replicate_groups[replicate_group]['peak_dir']+'align2rawsignal'
                    os.system('mkdir '+outdir)
                    n_reads_given=str(int(opts.subsample_reads)*1000000)
                    print 'here1--00000000000'
                    replicate_col1=re.sub('.tagAlign.gz','','Rep_0_'+replicate_group+'.mergedReplicates.tagAlign.gz')+'.subsampleTo'+n_reads_given+'.tagAlign.gz'
                    replicate_col2=re.sub('tagAlign.gz','wiggler',replicate_col1)
                    cur_pairfile=outdir+'/'+replicate_col2+'.pairFile'
                    os.system('echo "'+replicate_col1+'\t'+replicate_col2+'" > '+cur_pairfile)
                    cur_fragLen=outdir+'/'+replicate_col2+'.fragLen'
                    total_window=opts.total_window
                    print 'here 1.5'
                    print indir
                    cur_fragLen_value=open(indir+'/QC/'+replicate_col1+'SPP_table.txt','r').readlines()[0].strip().split()[2].split(',')[0]
                    print 'here2 -------------------'
                    print replicate_col1
                    os.system('echo "'+replicate_col1+'\t'+cur_fragLen_value+'\t'+total_window+'" > '+cur_fragLen)
                    ############################## for now hardcoded, should change this   
                    outdir_new=opts.align2rawsignal_bigwigdir
                    export_tmp_cmd='export TMP='+outdir_new+replicate_col2+'TMP'
                    align2rawsignal_cmd='/srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/wiggler/batchWiggler.sge.sh '+cur_pairfile+' '+cur_fragLen+' '+indir+' '+outdir_new+' TRUE 20 bg epanechnikov 5 0 '+seqdir+' '+umap+' '+opts.chrSizes
                    os.system(export_tmp_cmd)
                    print 'here3'
                    print align2rawsignal_cmd
                    os.system(align2rawsignal_cmd)
                else:
                    print 'Under construction to call on the not subsampled things'
            except:
                print "not running "+replicate_group


    if opts.step_to_perform=='mergePeaks':
        #Find merge groups
        if not opts.peaks_on_subsampled:
            sys.exit("Only supporting this for subsampled peaks at the moment")
        merge_groups={}
        n_reads_given=str(int(opts.subsample_reads)*1000000)
        for replicate_group in replicate_groups.keys():
            #We're merging narrowPeak files
            #Rep_0_NA19099_H3K4ME1.mergedReplicates.subsampleTo50000000.tagAlign_VS_INPUT_7YRI_GM19239_GM19238_GM18505_GM19240_GM19099_GM19193_GM18486.subsampleTo77000000.tagAlign_peaks.narrowPeak
            replicate_group_specific='Rep_0_'+replicate_group+'.mergedReplicates'+'.subsampleTo'+n_reads_given+'.tagAlign'+'_VS_'+os.path.basename(list(replicate_groups[replicate_group]['input_name'])[0]).strip('.gz')
            cur_merge_group=replicate_groups[replicate_group]['merge_group']
            if cur_merge_group not in merge_groups.keys():
                merge_groups[cur_merge_group]={}
                merge_groups[cur_merge_group]['files']=set()
                merge_groups[cur_merge_group]['merge_dir']=replicate_groups[replicate_group]['merge_dir']
            merge_groups[cur_merge_group]['files'].add(replicate_groups[replicate_group]['peak_dir']+replicate_group_specific+'_peaks.narrowPeak')
        for merge_group in merge_groups.keys():
            os.system('mkdir '+merge_groups[merge_group]['merge_dir'])
            merge_name=merge_groups[merge_group]['merge_dir']+'mergePeaks_'+merge_group+'.gz' #gzip
            #pooling_cmd='cat '+' '.join(merge_groups[merge_group]['files'])+' > '+merge_name+'_cat'
            #Run mergeBed - Sofia says: zcat $files | cut -f1-3 | sort -V | mergeBed -i stdin | gzip -c > $outfile
            load_module_cmd='module load bedtools/2.18.0'
            merging_cmd='zcat '+' '.join(merge_groups[merge_group]['files'])+' | cut -f1-3,8 | sort -V | bedtools merge -i stdin -nms | gzip > '+merge_name
            #Remove peaks that have bad pvalue in whole group. 5=-log10(min pvalue)
            maxMinPvalue_minusLog10='5'
            remove_weak_peaks='/srv/gs1/software/R/R-3.1.0/bin/Rscript /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/chromoVariationCode/keepPeaks_by_maxMinPvalue.R '+merge_name+' '+maxMinPvalue_minusLog10+' '+merge_name+'_Log10PvalueThreshold_'+maxMinPvalue_minusLog10+'.gz'
            #and also get rid of blacklisted
            blacklist_cmd=''
            ungzblacklist=''
            remove_weak_peaks_blacklist=''
            print 'merge name'
            print os.path.dirname(merge_name)
            if opts.blacklist!='':
                blacklist_site=os.path.dirname(merge_name)+'/'+os.path.basename(opts.blacklist)+'notgz'
                ungzblacklist='zcat '+opts.blacklist+' > '+blacklist_site
                blacklist_cmd='zcat '+merge_name+' | subtractBed -A -a stdin -b '+blacklist_site+' | gzip > '+merge_name+'removeBlacklist.gz'
                remove_weak_peaks_blacklist='/srv/gs1/software/R/R-3.1.0/bin/Rscript /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/chromoVariationCode/keepPeaks_by_maxMinPvalue.R '+merge_name+'removeBlacklist.gz'+' '+maxMinPvalue_minusLog10+' '+merge_name+'removeBlacklist.gz'+'_Log10PvalueThreshold_'+maxMinPvalue_minusLog10+'.gz'
            cmds=[load_module_cmd,merging_cmd,remove_weak_peaks,ungzblacklist,blacklist_cmd,remove_weak_peaks_blacklist]
            #print '\n'.join(cmds)
            qsub_a_command('qqqq'.join(cmds),merge_name+'_script.sh','qqqq',memory_number='3G')
        

def initialize_wiggler():
    cmds=[]
    cmds.append('export MCRROOT="/srv/gsfs0/projects/kundaje/users/akundaje/local/lib/mcr/2010b/v714"')
    cmds.append('export BACKUP_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}')
    cmds.append('LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/runtime/glnxa64')
    cmds.append('LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64')
    cmds.append('LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64')
    cmds.append('MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64')
    cmds.append('LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads')
    cmds.append('LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server')
    cmds.append('LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}')
    cmds.append('XAPPLRESDIR=${MCRROOT}/X11/app-defaults')
    cmds.append('export LD_LIBRARY_PATH')
    cmds.append('export XAPPLRESDIR')
    cmds.append('export MCR_CACHE_ROOT=${TMP}')
    cmds.append('export PATH=/srv/gsfs0/projects/kundaje/users/akundaje/projects/encode/code/signalgeneration/align2rawsignal/trunk/bin/:$PATH')
    return cmds

def given_samples_compute_mean_fragment_length(fragLen_files):
    fragment_length_column=2
    fragLen=0
    how_many=0
    for fragLen_file in fragLen_files:
        fragment_length=int(open(fragLen_file,'r').readlines()[0].split()[fragment_length_column].split(',')[0]) #take first value
        print fragLen_file
        print fragment_length
        how_many=how_many+1
        fragLen=fragLen+fragment_length
    mean_fragment_length=fragLen/how_many
    return mean_fragment_length

def qsub_a_command(cmd,shell_script_name,split_string=',',memory_number='20G'):
    f=open(shell_script_name,'w')
    cmds=cmd.split(split_string)
    for i in range(len(cmds)):
        #f.write("cmd"+str(i)+"='"+cmds[i]+"'"+'\n')                                                                                                                     
        #f.write('echo $cmd'+str(i)+'\n')                                                                                                                                
        #f.write('eval $cmd'+str(i)+'\n')
        f.write(cmds[i]+'\n') #just write the command
    f.close()
    #make runnable                                                                                                                                                      
    os.system('chmod 711 '+shell_script_name)
    #Qsub the script                                                                                                                                                    
    os.system("qsub -l mem_free="+memory_number+" -l h_vmem="+memory_number+" -l h_rt=20:00:00 -o "+shell_script_name+'.o'+' -e '+shell_script_name+'.e'+' '+shell_script_name)
        

main()
