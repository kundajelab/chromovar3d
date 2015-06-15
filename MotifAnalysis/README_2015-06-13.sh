#setup of data paths
setup_vars(){
  outmain=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-06-13/
  source_place=/srv/gsfs0/projects/kundaje/users/oursu/code/sequencingUtilities/peakCalling/chromoVariationCode/MotifAnalysis/MotifAnalysis_bashrc.bashrc
  source ${source_place}
  metadata=${outmain}/metadata.75unrelatedYRI
  TF2name=${outmain}/TFBS/ENCODE.hg19.TFBS.QC.metadata.jun2012_TFs_SPP_pooled.tsv
  olddir=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/
  CODE=/srv/gsfs0/projects/kundaje/users/oursu/code/chromovar3d/
}

#bring in some files from the previous analysis
mkdir -p ${outmain}/TFBS
cp ${olddir}/TFBS/TFpeaks_Gm12878.gz.named.gz ${outmain}/TFBS/
cp ${olddir}/metadata.75unrelatedYRI ${outmain}/

# SNPQTL files
#=============
# LOCAL
#=======
setup_vars ####
source_place=${CODE}/MotifAnalysis/MotifAnalysis_bashrc.bashrc
source ${source_place}
snpfile=/srv/gs1/projects/snyder/jzaugg/histoneQTL/hQTL_analysis/results/20140622_newPeaks/GWAS/hQTLs_
mkdir -p ${outmain}Local/SNPQTLmatrix

for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA
do
	echo ${hmark}
	cat ${snpfile}${hmark}.chr*.FDR0.1.txt | awk -F "\t" '{print $11"\t"$6"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$9"\t"$10"\t"$12}' |  head -n1 > ${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}
	cat ${snpfile}${hmark}.chr*.FDR0.1.txt | awk -F "\t" '{print $11"\t"$6"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$9"\t"$10"\t"$12}' | grep -v beta >> ${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}
	cat ${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark} | gzip > ${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}.gz
	rm ${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}
 done

# Make fasta files around the SNPs of interest
#=============================================
#Prepare SNPs
#============
source_place=${CODE}/MotifAnalysis/MotifAnalysis_bashrc.bashrc
source ${source_place}
motifLens=7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33
TFBSnamed=${outmain}/TFBS/TFpeaks_Gm12878.gz.named.gz
filters=${TFBSnamed}
metadata=${outmain}metadata.75unrelatedYRI

motifLens=19

for hmark in H3K4ME1 H3K27AC H3K4ME3 dhs RNA
do
    motifData=${outmain}/Motif_Data/MotifData_${hmark}
    mkdir -p ${motifData}
    SNPQTL=${outmain}Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}.gz
    OUT=${outmain}Local/Motif_Analysis/Motif_Analysis_${hmark}
    mkdir -p ${OUT}
    #Prepare SNPs and fasta files
    ${CODEDIR}/bin/MotifAnalysis_prepareSNPs.sh ${motifLens} ${filters} ${metadata} ${SNPQTL} ${motifData} ${source_place}
done



