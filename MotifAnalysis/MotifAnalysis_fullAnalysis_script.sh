#setup of data paths
setup_vars(){
  outmain=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-05-30/
  source_place=/srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/chromoVariationCode/MotifAnalysis/MotifAnalysis_bashrc.bashrc
  source ${source_place}
  metadata=${outmain}/metadata.75unrelatedYRI
  TF2name=${outmain}/TFBS/ENCODE.hg19.TFBS.QC.metadata.jun2012_TFs_SPP_pooled.tsv
  olddir=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/
}

###############################################
########## COPIED FROM ANALYSIS FROM 2015-01-13
###############################################
#commands done in the past already, and then just copied for this version of the analysis
#just see README_2015-01-13 for the commands done to make the data that I just copied to here

################################################################################################################
################################################################################################################

################################################
######### modified now 2015-05-30
################################################

# Get files
#==========
setup_vars ########
mkdir -p ${outmain}/Data/
mkdir -p ${outmain}/OverlapEnrichment/
mkdir -p ${outmain}/TFBS/
cp ${olddir}/Data/CombinedPeakAnno.gz ${outmain}/Data/CombinedPeakAnno.gz
cp -r ${olddir}/OverlapEnrichment/Peaks ${outmain}/OverlapEnrichment/Peaks
cp ${olddir}/TFBS/TFpeaks_Gm12878.gz ${outmain}/TFBS/TFpeaks_Gm12878.gz
cp ${olddir}/TFBS/ENCODE.hg19.TFBS.QC.metadata.jun2012_TFs_SPP_pooled.tsv ${TF2name}
#### no, i should remake this with the corrected PEER - cp ${olddir}/Data/CombinedPeakSignal.gz ${outmain}/Data/CombinedPeakSignal.gz
#run combinePeakSignal_2015-05-30.R
cp ${olddir}/TFBS/TFpeaks_Gm12878.gz.named.gz ${outmain}/TFBS/TFpeaks_Gm12878.gz.named.gz
# Metadata
setup_vars ########
initial_metadata=/srv/gs1/projects/snyder/jzaugg/histoneQTL/ChIPseq_alignment/data/metaData/fullMetaData_forAlignment.Unrelated75.keep
cat ${initial_metadata} | grep -v "#" | cut -f1,5 | sort | uniq | awk '{print $1"\t"$2"/"$1".maternal.fa\t"$2"/"$1".paternal.fa"}' > ${metadata} 

# SNPQTL files
#=============
# LOCAL
#=======
setup_vars ########
snpfile=/srv/gs1/projects/snyder/jzaugg/histoneQTL/hQTL_analysis/results/20140622_newPeaks/GWAS/hQTLs_
mkdir -p ${outmain}Local/SNPQTLmatrix
mkdir -p ${outmain}Local/SNPQTLmatrix/
#snpqtl files and tracks
for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA
do
  echo ${hmark}
  cat ${snpfile}${hmark}.chr*.FDR0.1.txt | awk '{print $11"\t"$6"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$9"\t"$10"\t"$12}' |  head -n1 > ${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}
  cat ${snpfile}${hmark}.chr*.FDR0.1.txt | awk '{print $11"\t"$6"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$9"\t"$10"\t"$12}' | grep -v beta >> ${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}
  cat ${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark} | gzip > ${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}.gz
  rm ${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}
  zcat ${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}.gz | grep -v "snp.position" | awk '{snp=$1+1}{print "chr"$2"\t"$1"\t"snp"\t"$2":"$1"-"snp"\t"$0}' | gzip > ${outmain}Local/SNPQTLmatrix/SNPQTLmatrixBed.${hmark}.gz
  pairwise_track ${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}.gz ${outmain}/OverlapEnrichment/Peaks/peaks.${hmark}.anno.gz
done
#compare the local snpQTL matrices with the ones I had before
# to see if i can directly use the scored snps
setup_vars ########
check_snpqtl_matrices(){
  cd ${outmain}/test/
  for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA;do echo ${hmark};  cp ../../2015-01-13/Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}.gz jan_SNPQTLmatrix.${hmark}.gz;cp ../../2015-05-30/Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}.gz may_SNPQTLmatrix.${hmark}.gz;done
  for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA;do echo ${hmark};  zcat jan_SNPQTLmatrix.${hmark}.gz | awk '{print $2":"$1"-"$4}' | sort | uniq > jan_snps_${hmark};zcat may_SNPQTLmatrix.${hmark}.gz | awk '{print $2":"$1"-"$4}' | sort | uniq > may_snps_${hmark};wc -l jan_snps_${hmark};wc -l may_snps_${hmark};comm -12 jan_snps_${hmark} may_snps_${hmark} | wc -l;done
}
check_snpqtl_matrices
# got as below. so ok to copy-paste the already scored fastas and motif scores
# over from 2015-01-13
#H3K4ME1
#3866725 jan_snps_H3K4ME1
#3866725 may_snps_H3K4ME1
#3866725
#H3K4ME3
#1261510 jan_snps_H3K4ME3
#1261510 may_snps_H3K4ME3
#1261510
#H3K27AC
#2272630 jan_snps_H3K27AC
#2272630 may_snps_H3K27AC
#2272630
#dhs
#3370762 jan_snps_dhs
#3370762 may_snps_dhs
#3370762
#RNA
#273522 jan_snps_RNA
#273522 may_snps_RNA
#273522

#Get fasta files and motif matches, and TF data
#===============================================
setup_vars ########
#copy fasta files and motif matches
cp -r ${olddir}/Motif_Data ${outmain}/
#copy the motif scores
motifdirs=$(ls ${olddir}/Motif_analysis_results/Motif_Analysis_RNA/)
for motifdir in ${motifdirs};
do
  echo $motifdir 
  for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA;
  do 
    mkdir -p ${outmain}/Motif_analysis_results/Motif_Analysis_${hmark}/${motifdir}/Scores
    cp ${olddir}/Motif_analysis_results/Motif_Analysis_${hmark}/${motifdir}/Scores/*Scores.compiled.gz ${outmain}/Motif_analysis_results/Motif_Analysis_${hmark}/${motifdir}/Scores/
  done
done 
#also bring over the TFBS peaks
mkdir -p ${outmain}/TFBS
cp -r ${olddir}/TFBS/unifiedPeaksByTF ${outmain}/TFBS/

#motif match the snpqtl matrices
#===============================
setup_vars ########
mkdir -p ${outmain}Local/SNPQTLmatrix/motifMatched/
filterSNPQTLmatrix_by_motifMatches (){
  snpqtl=$1
  snpqtlbed=$2
  motifm=$3
  out=$4
  source_place=$5

  mkdir -p $(dirname ${out})
  s=${out}.sh
  echo "zcat ${snpqtl} | head -n1 > ${out}" > ${s}
  echo "source ${source_place}" >> ${s}
  echo "bedtools intersect -u -a ${snpqtlbed} -b ${motifm} | sort | uniq | cut --complement -f1-4 >> ${out}" >> ${s}
  echo "rm ${out}.gz" >> ${s}
  echo "gzip ${out}" >> ${s}
  chmod 755 ${s}
  qsub ${s}
}
for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA
do
  echo ${hmark}
  motifs=$(ls ${outmain}/Motif_Data/motifMatches/*OverlapChIPseq.gz)
  for m in $motifs
  do
    echo $m
    filterSNPQTLmatrix_by_motifMatches ${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}.gz \
    ${outmain}Local/SNPQTLmatrix/SNPQTLmatrixBed.${hmark}.gz \
    $m ${outmain}Local/SNPQTLmatrix/motifMatched/SNPQTLmatrix.${hmark}.$(basename $m | sed 's/.gz//g') \
    ${source_place}
  done
done

#**********************************************************************

#Distal QTLs
#============
#Distal QTLs are here:
#/srv/gs1/projects/snyder/jzaugg/histoneQTL/hQTL_analysis/results/20140622_newPeaks/allQTLs_states.prox.frag.FDR10.allChr.2000kb.<MARK>.HiC_TH0.4.LDcorrected_0.2.states.rda
#need to process them for my use
setup_vars ########
mkdir -p ${outmain}/Distal/
#processed to be here
#/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-05-30/Distal/DistalQTLs.2015-05-30

#make the SNPQTL matrices
#bring the distal SNPQTL matrices
setup_vars ########
mkdir -p ${outmain}Distal/SNPQTLmatrix
distaldata=${outmain}/Distal/DistalQTLs.2015-05-30.ChrGeneDistal.PlusOverlappedRNA.txt
cat ${distaldata} | sort -k3b,3 > ${distaldata}.sorted
#make SNPQTL matrix by local peak
for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA
do
    #col2=distal, col3=local
    echo ${hmark}
    localsnpqtl=${outmain}/Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}.gz
    zcat ${localsnpqtl} | awk -v histmark=$hmark '{print histmark"_"$4"\t"$2"\t"$3"\t"$1}' | sort -k1b,1 | gzip > ${localsnpqtl}.sorted.gz
    echo "localpeak_chr_SNP_snp.position_chromo_gene" | tr "_" "\t" > ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.LocalPeakIs${hmark}
    zcat ${localsnpqtl}.sorted.gz | join -1 1 -2 3 - ${distaldata}.sorted | sed 's/ /\t/g' >> ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.LocalPeakIs${hmark}
    cat ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.LocalPeakIs${hmark} | gzip > ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.LocalPeakIs${hmark}.gz
    rm ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.LocalPeakIs${hmark}
    rm ${localsnpqtl}.sorted.gz
done
rm ${distaldata}.sorted

#make SNPQTL matrix by distal peak
cd ${outmain}Distal/SNPQTLmatrix/
zcat ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.LocalPeakIs*.gz | grep -v "snp.pos" | sed 's/[_]/\t/g' | sed 's/ /\t/g' | awk '{print $0>"SNPQTLmatrix.DistalPeakIs"$7}'
for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA
do
  echo "localpeak_chr_SNP_snp.position_chromo_gene" | tr "_" "\t" > ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.DistalPeakIs${hmark}_tmp
  cat ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.DistalPeakIs${hmark} | awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"_"$8}' >> ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.DistalPeakIs${hmark}_tmp
  cat ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.DistalPeakIs${hmark}_tmp | awk -v histmark=$hmark '{gsub(histmark"_","",$6)}1' | gzip > ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.DistalPeakIs${hmark}.gz
  rm ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.DistalPeakIs${hmark}_tmp
  rm ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.DistalPeakIs${hmark}

  zcat ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.DistalPeakIs${hmark}.gz | \
  sed 's/ /\t/g' | grep -v "snp.position" | awk '{snp=$4+1}{print "chr"$2"\t"$4"\t"snp"\t"$2":"$4"-"snp"\t"$0}' | \
  gzip > ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrixBed.DistalPeakIs${hmark}.gz 
  motifs=$(ls ${outmain}/Motif_Data/motifMatches/*OverlapChIPseq.gz)
  for m in $motifs
  do
    echo $m
    filterSNPQTLmatrix_by_motifMatches ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.DistalPeakIs${hmark}.gz \
    ${outmain}Distal/SNPQTLmatrix/SNPQTLmatrixBed.DistalPeakIs${hmark}.gz \
    $m ${outmain}Distal/SNPQTLmatrix/motifMatched/SNPQTLmatrix.DistalPeakIs${hmark}.$(basename $m | sed 's/.gz//g') \
    ${source_place}
  done
done
#**************************************************

#================================
# LOCAL + DISTAL ASSOCIATION TEST
#================================
associationTest (){
  motifName=$1
  localdistal=$2

  combSignal=${outmain}/Data/CombinedPeakSignal.gz
  source ${source_place}
  TFBSnamed=${outmain}/TFBS/TFpeaks_Gm12878.gz.named.gz
  filters=${TFBSnamed}
  metadata=${outmain}metadata.75unrelatedYRI
  for hmark in H3K4ME1 H3K27AC H3K4ME3 dhs RNA
  do
   mfile=${outmain}/Motif_Data/motifMatches/MotifMatch_$(echo ${motifName} | sed 's/motif.pouya.Motif.//g').bed.gz.OverlapChIPseq.gz
   #step2
   motifData=${outmain}/Motif_Data/MotifData_${hmark}
   #step3 
   OUT=${outmain}Motif_analysis_results/Motif_Analysis_${hmark}
   signalFile=/srv/gs1/projects/snyder/jzaugg/histoneQTL/hQTL_analysis/PeakData/hMat.norm.ALL.${hmark}.peer_lS_5.txt
   if [[ ${hmark} == "RNA" ]];
   then
    signalFile=/srv/gs1/projects/snyder/jzaugg/histoneQTL/hQTL_analysis/PeakData/hMat.norm.ALL.${hmark}.peer_lS_10.txt
   fi
   PWMdir=/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/data/motifPWM/matchesPerTF/correctedByPouyaControls/motifFreqMatrices/
   PWM=${PWMdir}${motifName}
   motifDir=/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/data/motifPWM/matchesPerTF/correctedByPouyaControls/motifFreqMatrices/matchesPerTF/
   if [[ ${localdistal} == "localanalysis" ]];
    then
      SNPQTL=${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}.gz
      SNPQTLmatched=${outmain}Local/SNPQTLmatrix/motifMatched/SNPQTLmatrix.${hmark}.MotifMatch_$(echo ${motifName} | sed 's/motif.pouya.Motif.//g').bed.OverlapChIPseq.gz
      ${CODEDIR}/bin/MotifAnalysis_associationTest.sh ${PWM} ${filters} ${metadata} ${SNPQTL} ${motifData} ${source_place} ${OUT} ${signalFile} ${mfile} "local" #${SNPQTLmatched}
    fi
   if [[ ${localdistal} == "distal" ]];
    then
      SNPQTL=${outmain}Distal/SNPQTLmatrix/SNPQTLmatrix.LocalPeakIs${hmark}.gz
      SNPQTLmatched=${outmain}Distal/SNPQTLmatrix/motifMatched/SNPQTLmatrix.DistalPeakIs${hmark}.MotifMatch_$(echo ${motifName} | sed 's/motif.pouya.Motif.//g').bed.OverlapChIPseq.gz
      ${CODEDIR}/bin/MotifAnalysis_associationTest.sh ${PWM} ${filters} ${metadata} ${SNPQTL} ${motifData} ${source_place} ${OUT} ${combSignal} ${mfile} distal #${SNPQTLmatched}
   fi
  done
 }
setup_vars ########
cd ${outmain}/Motif_Data/motifMatches/
motifs=$(ls MotifMatch_*.bed.gz.OverlapChIPseq.gz | sed 's/MotifMatch_//g' | sed 's/.bed.gz.OverlapChIPseq.gz//g')
for motif in $motifs 
do 
  #associationTest "motif.pouya.Motif."$motif localanalysis 
  associationTest "motif.pouya.Motif."$motif distal
done

################################################################################################################

####################################################################################
## MOTIF CORRELATION POSTPROCESSING and FIGURE
##########################

#Compile datasets across TFs
compileDataAcrossTFs (){
  hmark=$1

  mkdir -p ${outmain}/CrossTFanalysis

  #Put together locals
  #====================
  cat ${outmain}/Motif_analysis_results/Motif_Analysis_${hmark}/motif.pouya.Motif.*scanThresh0/Correlations/*motifMAF0.05.correlation.local.tab | head -n1 > ${outmain}/CrossTFanalysis/LocalAnalysis.motifMAF05.${hmark}
  #here we are taking all entries minus the column names
  cat ${outmain}/Motif_analysis_results/Motif_Analysis_${hmark}/motif.pouya.Motif.*scanThresh0/Correlations/*motifMAF0.05.correlation.local.tab | grep -v "spearman" >> ${outmain}/CrossTFanalysis/LocalAnalysis.motifMAF05.${hmark}
  cat ${outmain}/CrossTFanalysis/LocalAnalysis.motifMAF05.${hmark} | gzip > ${outmain}/CrossTFanalysis/LocalAnalysis.motifMAF05.${hmark}.gz
  rm ${outmain}/CrossTFanalysis/LocalAnalysis.motifMAF05.${hmark}
  #annotate with QTL type
  snpqtl_local=${outmain}/Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}.gz
  motif_res=${outmain}/CrossTFanalysis/LocalAnalysis.motifMAF05.${hmark}.gz
  zcat ${snpqtl_local} | awk '{print $3"_"$4"\t"$0}' | sort -k1b,1 > ${snpqtl_local}.sort
  zcat ${motif_res} | awk '{print $2"_"$3"\t"$0}' | sort -k1b,1 > ${motif_res}.sort
  echo "SNP_gene "$(zcat ${snpqtl_local} | head -n1)" "$(zcat ${motif_res} | head -n1) | sed 's/ /\t/g'  > ${motif_res}.annotated
  join -1 1 -2 1 ${snpqtl_local}.sort ${motif_res}.sort | sed 's/ /\t/g' >> ${motif_res}.annotated
  cat ${motif_res}.annotated | gzip > ${motif_res}.annotated.gz
  rm ${motif_res}.sort ${motif_res}.annotated 

  #Put together distals
  #====================
  cat ${outmain}/Motif_analysis_results/Motif_Analysis_${hmark}/motif.pouya.Motif.*scanThresh0/Correlations/*distal.tab | head -n1 > ${outmain}/CrossTFanalysis/DistalAnalysis.DistalPeakIs${hmark}
  #here we are taking all entries minus the column names
  cat ${outmain}/Motif_analysis_results/Motif_Analysis_${hmark}/motif.pouya.Motif.*scanThresh0/Correlations/*distal.tab | grep -v "spearman" >> ${outmain}/CrossTFanalysis/DistalAnalysis.DistalPeakIs${hmark}
  cat ${outmain}/CrossTFanalysis/DistalAnalysis.DistalPeakIs${hmark} | gzip > ${outmain}/CrossTFanalysis/DistalAnalysis.DistalPeakIs${hmark}.gz
  rm ${outmain}/CrossTFanalysis/DistalAnalysis.DistalPeakIs${hmark}
}
setup_vars ####
for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA
do
  echo "Compiling data across TFs for ${hmark}"
  compileDataAcrossTFs ${hmark} ${outmain}
done

#######################################################################################################################################################
# OVERLAP ENRICHMENT
#====================
setup_vars ########
DATA=${outmain}/OverlapEnrichment/
peaks=${DATA}/Peaks
enrichCode=/srv/gsfs0/projects/kundaje/users/oursu/code/genome_utils/Features/TFs/overlapEnrichment/overlapEnrichment.sh
AnnotateQTLpeaks () {
   snpqtl=$1
   zcat ${snpqtl} | awk '{if ($9=="pass") print $0}' | cut -f4 | sort | uniq | awk '{print $1"\tQTLpeak"}' > $(dirname ${snpqtl})/$(basename ${snpqtl} | sed 's/SNPQTLmatrix.//g' | sed 's/.gz//g').QTLpeaks
}
AnnotateDistalSourceQTLpeaks () {
   snpqtl=$1
   zcat ${snpqtl} | sed 's/ /\t/g' | sed 's/_/\t/g' | cut -f2 | sort | uniq | awk '{print $1"\tQTLpeak"}' > $(dirname ${snpqtl})/$(basename ${snpqtl} | sed 's/SNPQTLmatrix.//g' | sed 's/.gz//g').QTLpeaks_affectingDistalPeaks
}

# make regions files for peaks
# ${outmain}/OverlapEnrichment/Peaks/peaks.${hmark}.gz
for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA
do
    #make a bed file of peaks
    mkdir -p ${peaks}
    #make a file with QTL peaks
    snpqtl=${outmain}/Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}.gz
    AnnotateQTLpeaks ${snpqtl}
    #make a file with distalQTL peaks (these are the local QTLs that affect a distal peak)
    distal_snpqtl=${outmain}/Distal/SNPQTLmatrix/SNPQTLmatrix.LocalPeakIs${hmark}.gz
    AnnotateDistalSourceQTLpeaks ${distal_snpqtl}
done

overlapEnrichment_wrapper (){
  DATA=${outmain}/OverlapEnrichment/
  peaks=${DATA}/Peaks
  enrichCode=/srv/gsfs0/projects/kundaje/users/oursu/code/genome_utils/Features/TFs/overlapEnrichment/overlapEnrichment.sh
  #inputs
  local regions=$1
  local TFData=$2
  local TFData_regex=$3
  local outdir=$4
  local outpref=$5
  local window=$6
  local step=$7
  #we always aggregate over the peaks, and compute enrichment relative to the significant QTL peaks
  #either local or distal
  if [[ $step == "step1_computeOverlaps" ]];then
    for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA
    do
      #### TF overlap QTL peaks
      overall_command="${enrichCode} $(echo ${regions} | sed 's/HMARK/'${hmark}'/g') ${TFData} ${TFData_regex} ${outdir} ${window} ${peaks}/peaks.${hmark}.gz ${source_place}"
      ${overall_command} $(echo ${outpref} | sed 's/HMARK/'${hmark}'/g') computeOverlaps
    done
  fi

  if [[ $step == "step2_enrich" ]];then
    echo $(echo ${outpref} | sed 's/HMARK/'${hmark}'/g')
    echo "====================="
    #summarize overlaps into a matrix and compute overlap enrichment
    for hmark in RNA H3K4ME1 H3K4ME3 H3K27AC dhs
      do
      echo "Computing ${hmark}"
      #### TF overlap QTL peaks
      overall_command="${enrichCode} $(echo ${regions} | sed 's/HMARK/'${hmark}'/g') ${TFData} ${TFData_regex} ${outdir} ${window} ${peaks}/peaks.${hmark}.gz ${source_place} $(echo ${outpref} | sed 's/HMARK/'${hmark}'/g')"
      ${overall_command} overlapMatrix
      #enrichment analysis
      #local
      qtl_peaks=${outmain}/Local/SNPQTLmatrix/${hmark}.QTLpeaks
      ${overall_command} enrichmentAnalysis ${qtl_peaks}
      #distal
      qtl_peaks_affectingDistal=${outmain}/Distal/SNPQTLmatrix/LocalPeakIs${hmark}.QTLpeaks_affectingDistalPeaks
      ${overall_command} enrichmentAnalysis ${qtl_peaks_affectingDistal}
    done
  fi
}
#TFBS vs QTL peaks
cmd="overlapEnrichment_wrapper \
${peaks}/peaks.HMARK.gz \
${outmain}/TFBS/unifiedPeaksByTF/ \
MergedPeaks_ChIPseq_*.gz \
${DATA}/TFBS_overlap_QTLpeaks0kb \
TFBS_overlap_HMARKQTLpeaks0kb___ \
1"
${cmd} step1_computeOverlaps
${cmd} step2_enrich

# Motifs vs QTL peaks
cmd="overlapEnrichment_wrapper \
${peaks}/peaks.HMARK.gz \
${outmain}/Motif_Data/motifMatches/
MotifMatch_*scanThresh0.bed.gz.OverlapChIPseq.gz \
${DATA}/Motif_overlap_QTLpeaks0kb \
Motif_overlap_HMARKQTLpeaks0kb___ \
1"
${cmd} step1_computeOverlaps
${cmd} step2_enrich

# Disrupted motifs vs QTL peaks. The disrupted locals.
cmd="overlapEnrichment_wrapper \
${peaks}/peaks.HMARK.gz \
${outmain}/CrossTFanalysis/SNPsInMotifsBed/Local/CorrelatedMotifs/
correlatedMotifs.motif.pouya.Motif.*scanThresh0.gz  \
${DATA}/MotifCorrelatedLocal_overlap_QTLpeaks0kb \
MotifCorrelatedLocal_overlap_HMARKQTLpeaks0kb___ \
1"
${cmd} step1_computeOverlaps
${cmd} step2_enrich


###########################################################################################
# Extract signal for signal analysis
setup_vars ########
echo_matlab_things (){
  curScript=$1
  tempo=$2
  mkdir -p ${tempo}
  echo "source /srv/gsfs0/projects/kundaje/commonRepository/src/lab_bashrc" > ${curScript}
  echo "activateMCR" >> ${curScript}
  echo "export MCR_CACHE_ROOT=${tempo}" >> ${curScript}
  echo "export TMP=${tempo}" >> ${curScript}  
}

extractSignalByMark (){
  bedFile=$1
  outdir=$2
  mkdir -p ${outdir}
  mod=$3
  SignalFolder=$4
 
  cd $SignalFolder
  #hmarks
  indFiles=$(ls *mat)
  for indFile in $indFiles
  do
    echo $indFile
    echo $bedFile
    outFile=$outdir'MotifDisruptingPeakOf.'$mod'_'$indFile'.signal'
    curScript=$outFile'_extractSignal.sh'
    bedcopy=$outFile'_extractSignal'$mod'bed_uniqueSites.bed'
    if [ ! -e ${outFile} ];
      then
        zcat ${bedFile} | cut -f1-3 | sort | uniq | awk '{print "chr"$0}' > ${bedcopy}
        tempo=/srv/gsfs0/scratch/oursu/$(basename ${outFile})
        echo_matlab_things ${curScript} ${tempo}
        echo "extractSignal -i=${bedcopy} -t=${SignalFolder}/${indFile} -of=cagt -if=bed -o=${outFile} -sl=5000 -sr=5000" >> ${curScript}
        echo "cat ${outFile} | gzip > ${outFile}.gz" >> ${curScript}
        chmod 711 ${curScript}
        qsub -l h_vmem=50G -o ${curScript}.o -e ${curScript}.e ${curScript}
    fi
  done
  cd ${outdir}
}

#histone marks
hmark=H3K27AC
extractSignalByMark \
${outmain}/CrossTFanalysis/SNPsInMotifsBed/Local/CorrelatedMotifs/correlatedMotifs.${hmark}.allTFs.gz \
${outmain}/signal/ \
${hmark} \
/srv/gsfs0/projects/snyder/oursu/histoneQTL/ChIPseq_alignment/results/peakCalls/align2rawsignal/
#dnase
extractSignalByMark \
${outmain}/CrossTFanalysis/SNPsInMotifsBed/Local/CorrelatedMotifs/correlatedMotifs.${hmark}.allTFs.gz \
${outmain}/signal/DNase/ \
${hmark} \
/srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/subsampled/align2rawsignal/
#mnase
mnase_file=/srv/gsfs0/projects/kundaje/users/akundaje/projects/encode/data/byDataType/signal/jan2011/wigglerMatFiles/wgEncodeSydhGm12878NucleosomeRep0.bw30.norm5.rawsignal.mat
mkdir -p ${outmain}/signal/MNase
m_out=${outmain}/signal/MNase/mnase.signal
s=${m_out}.sh
echo_matlab_things ${s} /srv/gsfs0/scratch/oursu/
bedcopy=${m_out}bedcopy.bed
bedFile=${outmain}/CrossTFanalysis/SNPsInMotifsBed/Local/CorrelatedMotifs/correlatedMotifs.${hmark}.allTFs.gz
zcat ${bedFile} | cut -f1-3 | sort | uniq | awk '{print "chr"$0}' > ${bedcopy}
echo "extractSignal -i=${bedcopy} -t=${mnase_file} -of=cagt -if=bed -o=${m_out} -sl=5000 -sr=5000" >> ${s}
echo "cat ${m_out} | gzip > ${m_out}.gz" >> ${s}
echo "rm ${m_out}" >> ${s}
chmod 711 ${s}
qsub -l h_vmem=100G -o ${s}.o -e ${s}.e ${s}








