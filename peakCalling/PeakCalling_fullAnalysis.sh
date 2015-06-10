#### 1. Create metadata
#######################
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/src/construct_metaData_table_2014-02-19.py

#### 2. Split alignment metadata into 4 phases: SEQUENCED_Phase1, IMPUTED, chromoPaper, SEQUENCED_Phase3_CG
##############################################################################################################################
#1. SEQUENCED_Phase1
cat /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment | grep -w S > /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1
#2. IMPUTED
cat /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment | grep -w I > /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.IMPUTED
#3. FROM PREVIOUS PAPER
cat /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment | grep -w 'chromoPaper' > /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.chromoPaper
#4. SEQUENCED_Phase3_CG
cat /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment | grep -w 'seqCG' > /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase3_CG

#### 3. Create vcf files for each individual. 
#############################################
#Already done.

##### 4. Create personal genomes
################################
#1. SEQUENCED INDIVIDUALS
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/MAPPING_wrapper.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1 --step_to_perform createGenome
#2. IMPUTED INDIVIDUALS
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/MAPPING_wrapper.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.IMPUTED --step_to_perform createGenome
#3. FROM PREVIOUS PAPER
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/MAPPING_wrapper.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.chromoPaper --step_to_perform createGenome
#4. CG PHASE3
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/MAPPING_wrapper.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase3_CG --step_to_perform createGenome

####### 5. Make sure Fastq read names for PE files agree
#1. SEQUENCED INDIVIDUALS
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/check_read_names_agree_for_PE.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1 --out /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/checkFastqReadNames_for_PE_agreement.SEQUENCED_Phase1
#2. IMPUTED INDIVIDUALS
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/check_read_names_agree_for_PE.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.IMPUTED --out /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/checkFastqReadNames_for_PE_agreement.IMPUTED
#3. PREVIOUS PAPER
#No need to check, we will use the already mapped datasets.
#4. CG PHASE3
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/check_read_names_agree_for_PE.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase3_CG --out /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/checkFastqReadNames_for_PE_agreement.SEQUENCED_Phase3_CG

####### 6. Map reads
#1. SEQUENCED INDIVIDUALS
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/MAPPING_wrapper.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1 --step_to_perform align
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/MAPPING_wrapper.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1 --step_to_perform reconcile
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/MAPPING_wrapper.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1 --step_to_perform tagAlign
#2. IMPUTED
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/MAPPING_wrapper.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.IMPUTED --step_to_perform align
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/MAPPING_wrapper.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.IMPUTED --step_to_perform reconcile
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/MAPPING_wrapper.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.IMPUTED --step_to_perform tagAlign
#3. ChromoVar paper
#Already aligned. Just do tagAligns.
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/MAPPING_wrapper.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.chromoPaper --step_to_perform tagAlign
#4. CG PHASE3
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/MAPPING_wrapper.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase3_CG --step_to_perform align
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/MAPPING_wrapper.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase3_CG --step_to_perform reconcile
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/MAPPING_wrapper.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase3_CG --step_to_perform tagAlign

###### 7. Call peaks. Prepare metadata, run SPP, merge replicates, subsample, combine inputs and run MACS. Align2rawsignal.
#1. PREVIOUS PAPER
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/src/augment_metadata_for_peakCalls.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.chromoPaper
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.chromoPaper_for_peakCalling --step_to_perform SPP
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.chromoPaper_for_peakCalling --step_to_perform mergeReplicates
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.chromoPaper_for_peakCalling --step_to_perform subsample --subsample_reads 50
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.chromoPaper_for_peakCalling --subsample_reads 50 --step_to_perform MACS --peaks_on_subsampled
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.chromoPaper_for_peakCalling --subsample_reads 50 --peaks_on_subsampled --step_to_perform align2rawsignal

# === INPUT === Combine YRI inputs from previous study. 7 inputs, each contributes 11M reads. Concatenate them into 1 final input file. No need for SPP for input files.
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.chromoPaper_for_peakCalling --step_to_perform subsample --subsample_reads 11 --sample_names_to_do GM18486_INPUT,GM18505_INPUT,GM19099_INPUT,GM19193_INPUT,GM19238_INPUT,GM19239_INPUT,GM19240_INPUT
input_dir=/srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/Alignments/chromoVar_done_alignments/tagAlign/merged_replicates/subsampled
YRI_input=${input_dir}/INPUT_7YRI_GM19239_GM19238_GM18505_GM19240_GM19099_GM19193_GM18486.subsampleTo77000000.tagAlign.gz
zcat ${input_dir}/Rep_0_SNYDER_HG19_GM19239_INPUT.mergedReplicates.subsampleTo11000000.tagAlign.gz ${input_dir}/Rep_0_SNYDER_HG19_GM19238_INPUT.mergedReplicates.subsampleTo11000000.tagAlign.gz ${input_dir}/Rep_0_SNYDER_HG19_GM18505_INPUT.mergedReplicates.subsampleTo11000000.tagAlign.gz ${input_dir}/Rep_0_SNYDER_HG19_GM19240_INPUT.mergedReplicates.subsampleTo11000000.tagAlign.gz ${input_dir}/Rep_0_SNYDER_HG19_GM19099_INPUT.mergedReplicates.subsampleTo11000000.tagAlign.gz ${input_dir}/Rep_0_SNYDER_HG19_GM19193_INPUT.mergedReplicates.subsampleTo11000000.tagAlign.gz ${input_dir}/Rep_0_SNYDER_HG19_GM18486_INPUT.mergedReplicates.subsampleTo11000000.tagAlign.gz | gzip > ${YRI_input}
ln -s ${YRI_input} /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/Alignments/tagAlign/merged_replicates/subsampled/INPUT_7YRI_GM19239_GM19238_GM18505_GM19240_GM19099_GM19193_GM18486.subsampleTo77000000.tagAlign.gz

#2. SEQUENCED INDIVIDUALS
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/src/augment_metadata_for_peakCalls.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1_for_peakCalling --step_to_perform SPP
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1_for_peakCalling --step_to_perform mergeReplicates
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1_for_peakCalling --step_to_perform subsample --subsample_reads 50
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1_for_peakCalling --subsample_reads 50 --peaks_on_subsampled --step_to_perform align2rawsignal

#3. IMPUTED INDIVIDUALS
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/src/augment_metadata_for_peakCalls.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.IMPUTED
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.IMPUTED_for_peakCalling --step_to_perform SPP
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.IMPUTED_for_peakCalling --step_to_perform mergeReplicates
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.IMPUTED_for_peakCalling --step_to_perform subsample --subsample_reads 50
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.IMPUTED_for_peakCalling --subsample_reads 50 --step_to_perform MACS --peaks_on_subsampled 
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.IMPUTED_for_peakCalling --subsample_reads 50 --peaks_on_subsampled --step_to_perform align2rawsignal

#4. PHASE3 SEQ PEOPLE
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/src/augment_metadata_for_peakCalls.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase3_CG
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase3_CG_for_peakCalling --step_to_perform SPP
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase3_CG_for_peakCalling --step_to_perform mergeReplicates
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase3_CG_for_peakCalling --step_to_perform subsample --subsample_reads 50
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase3_CG_for_peakCalling --subsample_reads 50 --step_to_perform MACS --peaks_on_subsampled
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase3_CG_for_peakCalling --subsample_reads 50 --peaks_on_subsampled --step_to_perform align2rawsignal

#8. Merge all peak calling metadata to use for peak calling ---- HERE ---- /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.ALL_for_peakCalling
cat /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.*_for_peakCalling | grep '#' | head -n1 > /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.ALL_for_peakCalling
cat /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.*_for_peakCalling | grep -v '#' >> /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.ALL_for_peakCalling

#9. Merge ALL peaks.
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.ALL_for_peakCalling --subsample_reads 50 --peaks_on_subsampled --step_to_perform mergePeaks

#10. Extractsignal.
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.ALL_for_peakCalling --subsample_reads 50 --peaks_on_subsampled --step_to_perform extractSignal

#11. Merge extracted signal into 3 big (but zipped) matrices.
/srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/chromoVariationCode/mergeSignalValues_ALLindividuals.sh

# ------- PLOTS ---------
#------------------------
#1. SPP plot for all individuals
/srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/plots/2014-06-20/plotSPP_table.sh
#2. Merged peaks
/srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/chromoVariationCode/mergeSignalValues_ALLindividuals.sh

#========== Now, merge peaks using only the 75 individuals of interest, extract signal, make signal matrix
#get individuals
/home/oursu/devtools/R-3.0.2/bin/Rscript /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/src/select75individuals.R
#makes /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.Unrelated75.keep

#augment their metadata 
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/src/augment_metadata_for_peakCalls.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.Unrelated75.keep

#9.2 Merge peaks
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.Unrelated75.keep_for_peakCalling --subsample_reads 50 --peaks_on_subsampled --step_to_perform mergePeaks

#10.2 Extract signal
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/peakCalling.py --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.Unrelated75.keep_for_peakCalling --subsample_reads 50 --peaks_on_subsampled --step_to_perform extractSignal

#11.2 Merge extracted signal infto a matrix
#-------- HERE WE LEFT OFF
##modify this /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/chromoVariationCode/mergeSignalValues_ALLindividuals.sh







#==================== PAST
#PEAKS
#Need to divide metadata, so we can call peaks in batches
hmark=H3K4ME3
cat /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1_for_peakCalling | grep ${hmark} > /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1_for_peakCalling_${hmark}only
hmark=H3K4ME1
cat /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1_for_peakCalling | grep ${hmark} > /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1_for_peakCalling_${hmark}only
hmark=H3K27AC
cat /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1_for_peakCalling | grep ${hmark} > /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1_for_peakCalling_${hmark}only


##############

#2014-04-29
#Plot the CV, cluster peaks into variable and non using MAD
/srv/gs1/software/R/R-3.1.0/bin/Rscript plotCV_thresholdMAD.R /srv/gs1/projects/snyder/jzaugg/histoneQTL/peakAnnotation/peakLocations.H3K27ACvarScore.tab /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/plots/2014-04-29_H3K27AC H3K27AC /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/peakCalls/merged_peaks_H3K27AC/mergePeaks_H3K27AC.gzremoveBlacklist.gz_variablePeaks
/srv/gs1/software/R/R-3.1.0/bin/Rscript plotCV_thresholdMAD.R /srv/gs1/projects/snyder/jzaugg/histoneQTL/peakAnnotation/peakLocations.H3K4ME1varScore.tab /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/plots/2014-04-29_H3K4ME1 H3K4ME1 /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/peakCalls/merged_peaks_H3K4ME1/mergePeaks_H3K4ME1.gzremoveBlacklist.gz_variablePeaks
/srv/gs1/software/R/R-3.1.0/bin/Rscript plotCV_thresholdMAD.R /srv/gs1/projects/snyder/jzaugg/histoneQTL/peakAnnotation/peakLocations.H3K4ME3varScore.tab /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/plots/2014-04-29_H3K4ME3 H3K4ME3 /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/peakCalls/merged_peaks_H3K4ME3/mergePeaks_H3K4ME3.gzremoveBlacklist.gz_variablePeaks

#Combine this with definition of strong peaks for final set of peaks to be used for analysis
/srv/gs1/software/R/R-3.1.0/bin/Rscript getPeakIds.R /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/peakCalls/merged_peaks_H3K4ME1/mergePeaks_H3K4ME1.gzremoveBlacklist.gz_variablePeaks /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/peakCalls/merged_peaks_H3K4ME1/mergePeaks_H3K4ME1.gzremoveBlacklist.gz_Log10PvalueThreshold_5.gz
/srv/gs1/software/R/R-3.1.0/bin/Rscript getPeakIds.R /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/peakCalls/merged_peaks_H3K4ME3/mergePeaks_H3K4ME3.gzremoveBlacklist.gz_variablePeaks /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/peakCalls/merged_peaks_H3K4ME3/mergePeaks_H3K4ME3.gzremoveBlacklist.gz_Log10PvalueThreshold_5.gz 
/srv/gs1/software/R/R-3.1.0/bin/Rscript getPeakIds.R /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/peakCalls/merged_peaks_H3K27AC/mergePeaks_H3K27AC.gzremoveBlacklist.gz_variablePeaks /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/peakCalls/merged_peaks_H3K27AC/mergePeaks_H3K27AC.gzremoveBlacklist.gz_Log10PvalueThreshold_5.gz

/srv/gs1/software/R/R-3.1.0/bin/Rscript getPeakIds.R /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/peakCalls/merged_peaks_H3K4ME1/mergePeaks_H3K4ME1.gzremoveBlacklist.gz_variablePeaksNOT /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/peakCalls/merged_peaks_H3K4ME1/mergePeaks_H3K4ME1.gzremoveBlacklist.gz_Log10PvalueThreshold_5.gz
/srv/gs1/software/R/R-3.1.0/bin/Rscript getPeakIds.R /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/peakCalls/merged_peaks_H3K27AC/mergePeaks_H3K27AC.gzremoveBlacklist.gz_variablePeaksNOT /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/results/peakCalls/merged_peaks_H3K27AC/mergePeaks_H3K27AC.gzremoveBlacklist.gz_Log10PvalueThreshold_5.gz 
