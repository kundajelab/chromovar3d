setup_vars(){
source_file=/srv/gsfs0/projects/kundaje/users/oursu/code/genome_utils/Bedtools_like/Bedtools_for_pairs.sh
source ${source_file}
# datasets
#=========
motifdir=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-05-30
DATA=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-06-13
gwasdir=${DATA}/gwasOverlaps
qtldir=${DATA}/QTLs/
mkdir -p ${qtldir}
mkdir -p ${gwasdir}
#==================================================
# all marks and RNA annotations
annos=${DATA}/Annotations/
mkdir -p ${annos}
peakanno=${motifdir}/Data/CombinedPeakAnno.gz
# ENS annotation of TSSs (with gene symbols added if possible)
tss=${annos}/gencode.v19.annotation.PC.lincRNA.TSS.bed
regElements=${annos}/regElts.mergedHistonesDhs.bed
peaks2regElements=${annos}/peaks2regElts.mergedHistonesDhs.txt
peaks2TSS=${annos}/peaks2TSS.txt
TSSRegElt=${annos}/TSSRegEltNodeAnno.txt
code=/srv/gsfs0/projects/kundaje/users/oursu/code/chromovar3d/GWASnetworks
}

#bring in the GWAS files
cp -r /srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-05-30/GWAS_hits_Peyton ${DATA}/

#======================================
# Annotations
#======================================
make_annotations(){
  setup_vars
  #get a file of gencode TSSs. Need to have all genes from our dataset.
  python ${code}/gtf_to_geneTSSbed.py --gtf /srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/data/GENCODE_v19_2014-06-03/gencode.v19.annotation.PC.lincRNA.gtf --bed ${tss}
  #make the merged elements enhancers
  module load bedtools/2.21.0
  zcat -f ${peakanno} | grep -v RNA | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - | \
  awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > ${regElements}
  #map individual elements to the regulatory elements
  zcat -f ${peakanno} | grep -v RNA | cut -f1-4 | \
  bedtools window -w 0 -a - -b ${regElements} | cut -f4,8 | sort -k1b,1 > ${peaks2regElements}
  # map histonepeaks/dhs to the genes
  zcat -f ${peakanno} | cut -f1-4 | \
  bedtools window -w 0 -a - -b ${tss} | cut -f4,8,9 | sort -k1b,1 > ${peaks2TSS}
  #color nodes by TSS or reg elts
  zcat -f ${peaks2TSS} | awk '{print $2"\tTSS"}' > ${TSSRegElt}.tmp
  zcat -f ${regElements} | awk '{print $4"\tRegElt"}' >> ${TSSRegElt}.tmp
  echo "node_TSSRegElt" | sed 's/_/\t/g' > ${TSSRegElt}
  zcat -f ${TSSRegElt}.tmp | awk '!_[$1]++' >> ${TSSRegElt}
  rm ${TSSRegElt}.tmp
}
make_annotations
#=================

#==============
#make qtl data
#==============
setup_vars
zcat -f ${motifdir}/Local/SNPQTLmatrix/SNPQTLmatrix.*.gz | \
awk -F "\t" -v qtldirec=${qtldir} '{endsnp=$1}{beginsnp=$1-1}{if ($9=="pass") print $2"\t"beginsnp"\t"endsnp"\t"$3","$2":"beginsnp"-"endsnp",Local,"$10"_"$4",NA">qtldirec"/chr"$2".LocalQTLs.bed"}'
zcat -f ${motifdir}/Distal/DistalQTLs.2015-05-30 | \
awk -F "\t" -v qtldirec=${qtldir} '{beginsnp=$7-1}{endsnp=$7}{print "chr"$17"\t"beginsnp"\t"endsnp"\t"$1","$17":"beginsnp"-"endsnp",Distal,"$14","$2>qtldirec"/chr"$17".DistalQTLs.bed"}'
rm ${qtldir}/chrchr.DistalQTLs.bed

#also make a file with SNPs sorted by pvalue
#snp, peak, pvalue
#setup_vars
#zcat -f ${motifdir}/Local/SNPQTLmatrix/SNPQTLmatrix.*.gz | \
#awk -F "\t" '{if ($9=="pass") print $1"\t"$3"\t"$6"\t"$10"_"$4}' | grep -v "mod_gene" | \
#sort -g -k3 > ${qtldir}/pvalueSortedSNPspeakPeak.Local
#zcat -f ${motifdir}/Distal/DistalQTLs.2015-05-30 | \
#awk -F "\t" '{print $7"\t"$1"\t"$4"\t"$2}' | sort | uniq > ${qtldir}/pvalueSortedSNPspeakPeak.Distal
#zcat -f ${qtldir}/pvalueSortedSNPspeakPeak.Local ${qtldir}/pvalueSortedSNPspeakPeak.Distal | sort | uniq > ${qtldir}/pvalueSortedSNPspeakPeak
#python ${code}/SNPQTL_to_dict_bestSNPperPeak.py --snpqtl /srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-06-09/QTLs//pvalueSortedSNPspeakPeak --out /srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-06-09/QTLs//pvalueSortedSNPspeakPeak.dict
dict=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-06-09/QTLs//pvalueSortedSNPspeakPeak.dict
cp ${dict} ${DATA}/QTLs//pvalueSortedSNPspeakPeak.dict
#===============


#=======================================
# Expand in YRI LD the QTLs
#=======================================
expand_LD(){
  setup_vars ####
  local lQTL=$1
  local LDfile=/srv/gsfs0/projects/kundaje/users/pgreens/data/LD/YRI_geno.txt.gz
  local s=$(echo ${lQTL}expandLD.$(basename ${LDfile}).sh | sed 's/bed/gz/g')

  echo "module load tabix/0.2.6" > ${s}
  echo "python ${code}/QTLs_expandLD_fromTabix.py --LDfile ${LDfile} --infile ${lQTL}" >> ${s}
  chmod 755 ${s}
  qsub -o ${s}.o -e ${s} ${s}
}
setup_vars ####
for chromo in {1..22};
do 
  expand_LD ${qtldir}/chr${chromo}.LocalQTLs.bed
  expand_LD ${qtldir}/chr${chromo}.DistalQTLs.bed
done
#check them
cd ${qtldir}
for chromo in {1..22};do echo $chromo;zcat -f chr${chromo}.LocalQTLs.expandLD_YRI_geno.txt.new.bedbedbed | cut -f3 | sort | uniq | wc -l;zcat -f chr${chromo}.LocalQTLs.bed | cut -f3 | sort | uniq | wc -l;echo "======";done
for chromo in {1..22};do echo $chromo;zcat -f chr${chromo}.DistalQTLs.expandLD_YRI_geno.txt.new.bedbedbed | cut -f3 | sort | uniq | wc -l;zcat -f chr${chromo}.DistalQTLs.bed | cut -f3 | sort | uniq | wc -l;echo "======";done
#=========================================

#put them all together into 1 file
zcat -f ${qtldir}/chr*.DistalQTLs.expandLD_YRI_geno.txt.new.bedbedbed | \
awk '{print $9"\t"$10"\t"$11"\t"$4"\t"$13"\t"$12}' > ${qtldir}/ALL.DistalQTLs_expandYRILD.bed
zcat -f ${qtldir}/chr*.LocalQTLs.expandLD_YRI_geno.txt.new.bedbedbed | \
awk '{print $9"\t"$10"\t"$11"\t"$4"\t"$13"\t"$12}' > ${qtldir}/ALL.LocalQTLs_expandYRILD.bed
zcat -f ${qtldir}/ALL.DistalQTLs_expandYRILD.bed ${qtldir}/ALL.LocalQTLs_expandYRILD.bed > ${qtldir}/ALL.LocalandDistalQTLs_expandYRILD.bed
columns=chrLDSNP_startLDSNP_endLDSNP_QTL,chr:start-end,LocalDistal,LocalPeak,DistalPeak
echo "${columns}_R2_LDSNP" | \
sed 's/_/\t/g' > ${qtldir}/ALL.LocalandDistalQTLs_expandYRILD.bed.columns
#===========================

#==============
# Convert to our regulatory elements (merged histone + dhs), TSSs
#==========================================================
#convert to regulatory elements the name of the local peak and the name of the distal peak
LOCAL_PEAKCOL=7
DISTAL_PEAKCOL=8
columns=chrLDSNP_startLDSNP_endLDSNP_QTL,chr:start-end,LocalDistal,LocalPeak,DistalPeak
pref=${qtldir}/ALL.LocalandDistalQTLs_expandYRILD
zcat -f ${qtldir}/ALL.LocalandDistalQTLs_expandYRILD.bed | \
sed 's/,/\t/g' | sort -k${LOCAL_PEAKCOL}b,${LOCAL_PEAKCOL} | \
join -e "NA" -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 2.2 -a1  -1 ${LOCAL_PEAKCOL} -2 1 - ${peaks2regElements} | \
sed 's/ /\t/g' | \
sort -k${DISTAL_PEAKCOL}b,${DISTAL_PEAKCOL} | \
join -e "NA" -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 2.2 -a1  -1 ${DISTAL_PEAKCOL} -2 1 - ${peaks2regElements} | \
sed 's/ /\t/g' > ${pref}.2regElts.bed 
echo "${columns},R2,LDSNP,RegEltLocal,RegEltDistal" | \
sed 's/_/\t/g' | sed 's/,/\t/g' > ${pref}.2regElts.bed.columns
#convert to TSS the name of the local peak and the name of the distal peak
zcat -f ${pref}.2regElts.bed | \
sed 's/,/\t/g' | sort -k${LOCAL_PEAKCOL}b,${LOCAL_PEAKCOL} | \
join -e "NA" -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 2.2 -a1 -1 ${LOCAL_PEAKCOL} -2 1 - ${peaks2TSS} | \
sed 's/ /\t/g' | \
sort -k${DISTAL_PEAKCOL}b,${DISTAL_PEAKCOL} | \
join -e "NA" -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 2.2 -a1 -1 ${DISTAL_PEAKCOL} -2 1 - ${peaks2TSS} | \
sed 's/ /\t/g' > ${pref}.2regElts.2TSS.bed 
echo "${columns},R2,LDSNP,RegEltLocal,RegEltDistal,TSSLocal,TSSDistal" | \
sed 's/_/\t/g' | sed 's/,/\t/g' > ${pref}.2regElts.bed.2TSS.columns

#====================================
# Overlap with GWAS and make networks
#====================================
make_detailed_net(){
  local infile=$1
  local outfile=$2
  local pval_thresh=$3
  local LDthresh=$4
  local script=$5

  passFilters=${outfile}.LD${LDthresh}.p${pval_thresh}.txt

  echo "zcat -f ${infile} | \
  awk '{if ((\$19<=${pval_thresh}) && (\$9>=${LDthresh})) print \$0}' > ${passFilters}" >> ${script}
  echo "python ${code}/network_summary2.py --infile ${passFilters} --p ${pval_thresh}" >> ${script}
  echo "zcat -f ${passFilters}.net.pval${pval_thresh} | sort | uniq > ${passFilters}.net.pval${pval_thresh}.uniq" >> ${script}
}

setup_vars ####
columns=chrLDSNP_startLDSNP_endLDSNP_QTL,chr:start-end,LocalDistal,LocalPeak,DistalPeak
GWASPeyton=${DATA}/GWAS_hits_Peyton/
gwases=$(echo $(ls ${GWASPeyton}/expanded_LD_geno/rsq_0.8/*/*bed))
pref=${qtldir}/ALL.LocalandDistalQTLs_expandYRILD
bed_interest=${pref}.2regElts.2TSS.bed
#gwas=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-05-30/GWAS_hits_Peyton//expanded_LD_geno/rsq_0.8/roadmap/EUR.CD_pruned_rsq_0.8_expanded_rsq_0.8.bed

for gwas in ${gwases};
do
 out=${gwasdir}/overlapQTL_$(basename ${gwas} | sed 's/.bed//g')/overlapQTL_$(basename ${gwas}).LD08_p5.bed
 mkdir -p $(dirname ${out})
 s=${out}.sh
 GWAS_SNP_COL=18
 #overlap snps with the gwas
 echo "module load bedtools/2.21.0" > ${s}
 echo "source ${source_file}" >> ${s}
 #echo "zcat -f $(echo ${gwas} | sed 's/no_TAG_SNPs_in_LD_w_SNPs_in_LD.bed/no_TAG_SNPs_in_LD.bed/g') | sort -k4b,4 > ${out}.pvals" >> ${s}
 echo "cat $(echo ${gwas} | sed 's/expanded_LD_geno/pruned_LD_geno/g' | sed 's/_expanded_rsq_0.8//g') | sort -k4b,4 > ${out}.pvals" >> ${s}
 echo "bedtools window -w 0 -a ${bed_interest} -b ${gwas} | \
 sort -k${GWAS_SNP_COL}b,${GWAS_SNP_COL} | \
 join -1 4 -2 ${GWAS_SNP_COL} -o 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 2.10 2.11 2.12 2.13 2.14 2.15 2.16 2.17 2.18 1.5  ${out}.pvals - | \
 sed 's/ /\t/g' > ${out}.overlapGWAS.bed" >> ${s} #this gives us the QTL SNPs to keep (we'll keep them + all their LD)
 echo "echo ${columns},R2,LDSNP,RegEltLocal,RegEltDistal,TSSLocal,TSSDistal_GWASchr_GWASstart_GWASend_GWAStagSNPID_GWASassociationPval | \
 sed 's/_/\t/g' | sed 's/_/\t/g' > ${out}.overlapGWAS.bed.columns" >> ${s}
 make_detailed_net ${out}.overlapGWAS.bed ${out}.overlapGWAS 0.00001 0.8 ${s}
 chmod 755 ${s}
 qsub -o ${s}.o -e ${s}.e ${s}
done


### TABLE
#######################################
LD=0.8
p=0.00001

setup_vars ###
GWASPeyton=${DATA}/GWAS_hits_Peyton/
gwases=$(echo $(ls ${GWASPeyton}/expanded_LD_geno/rsq_0.8/*/*bed))
for gwas in ${gwases};
do
  echo $gwas 
  out=${gwasdir}/overlapQTL_$(basename ${gwas} | sed 's/.bed//g')/overlapQTL_$(basename ${gwas}).LD08_p5.bed
  zcat -f ${out}.overlapGWAS.LD${LD}.p${p}.txt | \
  awk '{print $18"SPLIT"$7"SPLIT"$8"SPLIT"$11"SPLIT"$12"SPLIT"$13"SPLIT"$14"SPLIT"$6"\t"$4}' | \
  sort -k1,1 | bedtools groupby -g 1 -c 2 -o distinct -i stdin | \
  sed 's/SPLIT/\t/g' > ${out}.overlapGWAS.LD${LD}.p${p}.GWAStable.txt_init
  zcat -f ${out}.overlapGWAS.LD${LD}.p${p}.GWAStable.txt_init | awk -v gwasname=$(echo $(basename ${gwas} | sed 's/_pruned_rsq_0.8_expanded_rsq_0.8.bed//g')) '{print $0"\t"gwasname}' > ${out}.overlapGWAS.LD${LD}.p${p}.GWAStable.txt
  rm ${out}.overlapGWAS.LD${LD}.p${p}.GWAStable.txt_init
done

table_file=${gwasdir}/GWAStable.pval5.simpleTable.txt
echo "GWAStagSNP_LocalPeak_DistalPeak_LocalRegElt_DistalRegElt_LocalGene_DistalGene_QTLType_QTLSNPs_Disease" | \
sed 's/_/\t/g' > ${table_file}
zcat -f ${gwasdir}/overlapQTL_*/*overlapGWAS.LD0.8.p0.00001.GWAStable.txt >> ${table_file}












