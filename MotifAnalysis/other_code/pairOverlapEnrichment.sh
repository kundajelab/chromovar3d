OVERLAP_CODE=/srv/gsfs0/projects/kundaje/users/oursu/code/genome_utils/Features/TFs/overlapEnrichment/
export PATH=$PATH:$OVERLAP_CODE

# First, do the simple TF overlap and create element-TF matrix
regionFile=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Data/CombinedPeakAnno.gz
TFData=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/TFBS/unifiedPeaksByTF/
TFData_fileRegex=MergedPeaks_ChIPseq_*.gz
outdir=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/pairAnalysis/TF_TF/
overlapWindow=0
aggregatedOver=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits/Annotations/regElts.mergedHistonesDhs.bed
bashrc=/srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/chromoVariationCode/MotifAnalysis/MotifAnalysis_bashrc.bashrc
out_prefix=OverlapPeaks_aggregateOver_regElts.mergedHistonesDhs.bed

#interactions=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/pairAnalysis/data/GWAS.DistalQTL.RegElts.bedpe
#get the interactions bedpe file, which will be a collection of the distal QTLs. 
#Can use directly the QTL peaks, they will automatically be aggregated across the reg elements
#
#distals=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits/QTLs/chr*.DistalQTLs.bed
#zcat -f ${distals} 


overlapEnrichment.sh ${regionFile} ${TFData} ${TFData_fileRegex} ${outdir} ${overlapWindow} ${aggregatedOver} ${bashrc} ${out_prefix} "computeOverlaps"
overlapEnrichment.sh ${regionFile} ${TFData} ${TFData_fileRegex} ${outdir} ${overlapWindow} ${aggregatedOver} ${bashrc} ${out_prefix} "overlapMatrix"

#let's get the interactions we have in a bedpe file
echo "chr_start_end_chr_start_end" | sed 's/_/\t/g' > ${interactions}
zcat -f /srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits/gwasOverlaps/GWAStable.total.txt | \
grep Distal | cut -f11,12 | sed 's/:/\t/g' | sed 's/-/\t/g' | grep -v RegElt | grep -v NA | sort | uniq >> ${interactions}

# Now, take the interactions, and add columns corresponding to the TFs
overlapEnrichment.sh ${regionFile} ${TFData} ${TFData_fileRegex} ${outdir} ${overlapWindow} ${aggregatedOver} ${bashrc} ${out_prefix} "mergeOverlapMatrixByBedpe" NA ${interactions}
overlapEnrichment.sh ${regionFile} ${TFData} ${TFData_fileRegex} ${outdir} ${overlapWindow} ${aggregatedOver} ${bashrc} ${out_prefix} "enrichmentAnalysisByBedpe" NA ${interactions}






#now, pick the columns corresponding to the 2 TFs, and count the occurrence of the pairs
tf1_1=$(zcat -f ${outdir}/${out_prefix}.site1site2Data.gz  | head -n1 | sed 's/\t/\n/g' | grep -n site1${tf1} | cut -f1 -d ":")
tf1_2=$(zcat -f ${outdir}/${out_prefix}.site1site2Data.gz  | head -n1 | sed 's/\t/\n/g' | grep -n site2${tf1} | cut -f1 -d ":")
tf2_1=$(zcat -f ${outdir}/${out_prefix}.site1site2Data.gz  | head -n1 | sed 's/\t/\n/g' | grep -n site1${tf2} | cut -f1 -d ":")
tf2_2=$(zcat -f ${outdir}/${out_prefix}.site1site2Data.gz  | head -n1 | sed 's/\t/\n/g' | grep -n site2${tf2} | cut -f1 -d ":")
zcat -f ${outdir}/${out_prefix}.site1site2Data.gz | awk -v tf11=${tf1_1} -v tf12=${tf1_2} -v tf21=${tf2_1} -v tf22=${tf2_2} -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$tf11"\t"$tf12"\t"$tf21"\t"$tf22}' | \


> ${outdir}/${out_prefix}.${tf1}.${tf2}

enrichmentAnalysisByBedpe

zcat -f ${outdir}/${out_prefix}.site1site2Data.gz | awk -v tf11=${tf1_1} -v tf12=${tf1_2} -v tf21=${tf2_1} -v tf22=${tf2_2} -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$tf11"\t"$tf12"\t"$tf21"\t"$tf22}' | awk '{if (($7>0 && $10>0) || ($8>0 && $9>0)) print $0"\tpair";else print $0"\tNotPair"}'