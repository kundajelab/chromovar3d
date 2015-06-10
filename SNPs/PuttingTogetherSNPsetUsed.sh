
#==================================================
#Get VCF data for our 109 individuals in the study.
#==================================================

#Subset individuals from Phase1. No chrY.
#========================================
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/src/subset_people_from_vcf.py --out /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/mergedVCF/phase1_subset. --out_vcf /srv/gsfs0/projects/kundaje/users/oursu/histoneQTLproject/out_vcf --chromo 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X --vcf_start /srv/gs1/projects/kundaje/oursu/Alignment/data/1000Genomes/VCF/ALL.chr --vcf_end .phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase1

#Subset sequenced individuals from Phase3. No chrY.
#==================================================
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/src/subset_people_from_vcf.py --out /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/mergedVCF/phase3seqCG_subset. --out_vcf /srv/gsfs0/projects/kundaje/users/oursu/histoneQTLproject/out_vcf --vcf_start /srv/gs1/projects/kundaje/oursu/Alignment/data/1000Genomes/VCF/ALL.chr --vcf_end .cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.genotypes.vcf.gz --chromo 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase3_CG

#Merge imputed people. No chrY. chrX has only phasing data, no imputation.
#===========================================
#Get imputation data to vcf
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/imputation_to_vcf/imputation_to_vcf.py --do_phasing 
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/imputation_to_vcf/imputation_to_vcf.py --do_score
python /srv/gsfs0/projects/kundaje/users/oursu/code/personalGenomeAlignment/imputation_to_vcf/imputation_to_vcf.py --do_imputed
#VCF: Imputed people
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/src/merge_people_into_vcf.py --out /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/mergedVCF/mergedImputed. --vcf_start /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/imputed/VCF/ --vcf_end .IMPUTED.vcf --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.IMPUTED
#VCF: Imputed people also sequenced with CG
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/src/merge_people_into_vcf.py --out /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/mergedVCF/mergedIMPUTED.phase3seqCG --metadata /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.SEQUENCED_Phase3_CG --vcf_start /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/imputed/VCF/ --vcf_end .IMPUTED.vcf

#Merge vcfs from individuals from previous paper
================================================
#1. Create vcfs based on previous paper, for SNPs sequenced in Phase1.
#Make 1 VCF with all locations in Phase1, which will be used for querying personal genomes from chromoPaper.
#cat /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/mergedVCF/phase1_subset.chr*.vcf | grep -v '#' | sed 's/chr//g' | cut -f1-10 > /srv/gs1/projec#ts/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/VCF/Phase1_Positions.vcf
#VCF files for chromoPaper

#2. Split VCFs into chromosomes
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/chromoVariationCode/merge_chromoPaper_VCFs.py --vcf_start /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/VCF/chromoPaper --vcf_end .vcf --people NA18486,NA19193,NA19238,NA19239,NA19240 --out /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/chromoPaper/VCF --step splitChromo
#3. Merge each chromosome into 1 file.
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/chromoVariationCode/merge_chromoPaper_VCFs.py --vcf_start /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/VCF/chromoPaper --vcf_end .vcf --people NA18486,NA19193,NA19238,NA19239,NA19240 --step mergeVCF

#Merge all 4 large VCFs into 1. Clean up names, and replace imputation results with sequencing, if these are available.
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/chromoVariationCode/merge_all_VCFs_cleanNames.py --indir /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/mergedVCF/ --phase1 phase1_subset.chr --phase3_imp mergedIMPUTED.phase3seqCGchr --phase3_cg phase3seqCG_subset.chr --imp mergedImputed.chr --chromoPaper chromoPaper.chr --out /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/mergedVCF/all_YRI_compiled/all_YRI_compiled --step merge --chromo 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X
#Clean up names and prioritize sequencing for individuals who are sequenced in Phase3.
python /srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/chromoVariationCode/merge_all_VCFs_cleanNames.py --indir /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/mergedVCF/ --phase1 phase1_subset.chr --phase3_imp mergedIMPUTED.phase3seqCGchr --phase3_cg phase3seqCG_subset.chr --imp mergedImputed.chr --chromoPaper chromoPaper.chr --out /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/mergedVCF/all_YRI_compiled/all_YRI_compiled --step clean --chromo 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X

#One large VCF with all VCF data
#================================
#This is for making personal genomes
cat /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/mergedVCF/all_YRI_compiled/all_YRI_compiled.chr*.mergedAcrossSNPSources.final.vcf > /srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/genomes/mergedVCF/all_YRI_compiled/all_YRI_compiled.ALL_CHROMOSOMES.mergedAcrossSNPSources.final.vcf