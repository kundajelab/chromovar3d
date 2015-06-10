setup_vars(){
source_file=/srv/gsfs0/projects/kundaje/users/oursu/code/genome_utils/Bedtools_like/Bedtools_for_pairs.sh
source ${source_file}
# datasets
#=========
motifdir=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-05-30
DATA=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-05-30/
gwasdir=${DATA}/gwasOverlaps
qtldir=${DATA}/QTLs/
mkdir -p ${qtldir}
mkdir -p ${gwasdir}
distal_qtls=${motifdir}/Distal/DistalQTLs.2015-05-30
local_qtls=${motifdir}/Local/LocalQTLs.2015-02-25
# all marks and RNA annotations
annos=${DATA}/Annotations/
mkdir -p ${annos}
peakanno=${motifdir}/Data/CombinedPeakAnno.gz
# ENS annotation of TSSs (with gene symbols added if possible)
peakanno_RNA=/srv/gs1/projects/snyder/jzaugg/histoneQTL/peakAnnotation/peakLocations.ALL.RNA.txt
tss=${annos}/pointTSS.bed
ens2symbols=/srv/gsfs0/projects/snyder/oursu/histoneQTL/networks/results/ensgenenames2.sorted
regElements=${annos}/regElts.mergedHistonesDhs.bed
peaks2regElements=${annos}/peaks2regElts.mergedHistonesDhs.txt
peaks2TSS=${annos}/peaks2TSS.txt
TSSRegElt=${annos}/TSSRegEltNodeAnno.txt
}


#==========================
#Annotations
#==========================
#make TSS annotations
setup_vars ####
zcat -f ${peakanno_RNA} | sed 's/ /\t/g' | awk '{start=$3+1500}{end=$3+1501}{print $2"\t"start"\t"end"\t"$5}' | \
sort -k4b,4 | join -1 4 -2 1 -a1 -o 1.1 1.2 1.3 0 2.2 - ${ens2symbols} | grep -v gene_id | sed 's/ /\t/g' > ${tss}
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
zcat -f ${peaks2TSS} | awk '{print $3"\tTSS"}' > ${TSSRegElt}.tmp
zcat -f ${regElements} | awk '{print $4"\tRegElt"}' >> ${TSSRegElt}.tmp
echo "node_TSSRegElt" | sed 's/_/\t/g' > ${TSSRegElt}
zcat -f ${TSSRegElt}.tmp | awk '!_[$1]++' >> ${TSSRegElt}
rm ${TSSRegElt}.tmp
#===========================

#==============
#make qtl data
#==============
zcat -f ${motifdir}/Local/SNPQTLmatrix/SNPQTLmatrix.*.gz | \
awk -F "\t" -v qtldirec=${qtldir} '{endsnp=$1+1}{if ($9=="pass") print $2"\t"$1"\t"endsnp"\t"$3","$2":"$1"-"endsnp",Local,"$10"_"$4",NA">qtldirec"/chr"$2".LocalQTLs.bed"}'
zcat -f ${motifdir}/Distal/DistalQTLs.2015-05-30 | \
awk -F "\t" -v qtldirec=${qtldir} '{plusone=$7+1}{print "chr"$17"\t"$7"\t"plusone"\t"$1","$17":"$7"-"plusone",Distal,"$14","$2>qtldirec"/chr"$17".DistalQTLs.bed"}'
rm ${qtldir}/chrchr.DistalQTLs.bed
#===============

#======================
# LD expansion of QTLs
#====================== 
#Make a bed file with both types of QTLs, and expand with LD
columns=chrLDSNP_startLDSNP_endLDSNP_QTL,chr:start-end,LocalDistal,LocalPeak,DistalPeak
for chromo in {1..22};
do
 #Distals
 out=${qtldir}/chr${chromo}.DistalQTLs_expandYRILD.bed
 s=${out}.sh
 echo "python /srv/gsfs0/projects/kundaje/users/oursu/code/ChromatinVariationCode/gwas/expandLD_by_column.py --out ${out}.tmp --chromo ${chromo} --snpcol 4 --input ${qtldir}/chr${chromo}.DistalQTLs.bed" > ${s}
 echo "echo \"${columns}_R2_LDSNP\" | sed 's/_/\t/g' > ${out}.columns" >> ${s}
 echo "zcat -f ${out}.tmp | cut --complement -f4-6 > ${out}" >> ${s}
 echo "rm ${out}.tmp" >> ${s}
 chmod 711 ${s}
 qsub -o ${s}.o -e ${s}.e ${s}
 #Locals
 out=${qtldir}/chr${chromo}.LocalQTLs_expandYRILD.bed
 s=${out}.sh
 echo "python /srv/gsfs0/projects/kundaje/users/oursu/code/ChromatinVariationCode/gwas/expandLD_by_column.py --out ${out}.tmp --chromo ${chromo} --snpcol 4 --input ${qtldir}/chr${chromo}.LocalQTLs.bed" > ${s}
 echo "echo \"${columns}_R2_LDSNP\" | sed 's/_/\t/g' > ${out}.columns" >> ${s}
 echo "zcat -f ${out}.tmp | cut --complement -f4-6 > ${out}" >> ${s}
 echo "rm ${out}.tmp" >> ${s}
 chmod 711 ${s}
 #qsub -o ${s}.o -e ${s}.e ${s}
done

#=========================
# Alternate LD expansion
#==========================
module load tabix/0.2.6
LDdir=
#*************************************
#*************************************

#put them all together into 1 file
#zcat -f ${qtldir}/chr*.DistalQTLs_expandYRILD.bed | grep rs*  > ${qtldir}/ALL.DistalQTLs_expandYRILD.bed
#zcat -f ${qtldir}/chr*.LocalQTLs_expandYRILD.bed | grep rs* > ${qtldir}/ALL.LocalQTLs_expandYRILD.bed
zcat -f ${qtldir}/chr*.DistalQTLs_expandYRILD.bed > ${qtldir}/ALL.DistalQTLs_expandYRILD.bed
zcat -f ${qtldir}/chr*.LocalQTLs_expandYRILD.bed > ${qtldir}/ALL.LocalQTLs_expandYRILD.bed

zcat -f ${qtldir}/ALL.DistalQTLs_expandYRILD.bed ${qtldir}/ALL.LocalQTLs_expandYRILD.bed > ${qtldir}/ALL.LocalandDistalQTLs_expandYRILD.bed
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

#========================
# GWAS
#=========================
make_net (){
  local infile=$1
  local outfile=$2
  local pval_thresh=$3
  local LDthresh=$4
  local script=$5

  annos=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits/Annotations/
  peaks2TSS=${annos}/peaks2TSS.txt
  tss=${annos}/pointTSS.bed
  regElements=${annos}/regElts.mergedHistonesDhs.bed

  echo "zcat -f ${infile} | grep Local | \
  awk '{if ((\$19<=${pval_thresh}) && (\$9>=${LDthresh})) print \$0}' | \
  awk '{n1=\$18;n2=\$11}{if (\$13!=\"NA\") n2=\$13}{print n1\"\t\"n2\"\t\"\$6\"\t\"\$7}' | \
  sed 's/_/\t/g' | \
  sed 's/H3K4ME1/chromatin/g' | sed 's/H3K4ME3/chromatin/g' | sed 's/H3K27AC/chromatin/g' | sed 's/dhs/chromatin/g' | \
  cut -f1-4 | sort | uniq > ${outfile}" >> ${script}

  echo "zcat -f ${infile} | grep Distal | \
  awk '{if ((\$19<=${pval_thresh}) && (\$9>=${LDthresh})) print \$0}' | \
  awk '{n1=\$11;n2=\$12}{if (\$13!=\"NA\") n1=\$13}{if (\$14!=\"NA\") n2=\$14}{print n1\"\t\"n2\"\t\"\$6\"\t\"\$8}' | \
  sed 's/_/\t/g' | \
  sed 's/H3K4ME1/chromatin/g' | sed 's/H3K4ME3/chromatin/g' | sed 's/H3K27AC/chromatin/g' | sed 's/dhs/chromatin/g' | \
  cut -f1-4 | sort | uniq >> ${outfile}" >> ${script}

  #annotate the genes with their symbol
  echo "zcat -f ${peaks2TSS} | sort -k2b,2 > ${outfile}.TSSsorted" >> ${script}
  echo "zcat -f ${outfile} | \
  sort -k1b,1 | join -1 2 -2 1 -a2 -o 2.2 2.3 2.4 1.3 2.1 ${outfile}.TSSsorted - | sed 's/ /\t/g' | awk '{print \$4\"\t\"\$1\"\t\"\$2\"\t\"\$3}' | \
  sort -k2b,2 | join -1 2 -2 2 -a2 -o 2.1 2.3 2.4 1.3 2.2 ${outfile}.TSSsorted - | sed 's/ /\t/g' | awk '{print \$1\"\t\"\$4\"\t\"\$2\"\t\"\$3}' | \
  sort | uniq > ${outfile}.withGeneNames.tmp" >> ${script}
  #annotate the genes with coordinates, so that we can compute distances between elements in our networks
  echo "zcat -f ${tss} | cut -f1-3,5 | zcat -f - ${regElements} | cut -f1-4 | sort -k4b,4 > ${outfile}.allElementsBed" >> ${script}
  echo "zcat -f ${outfile}.withGeneNames.tmp | grep Distal | sort -k1b,1 | \
  join -1 1 -2 4 -a1 -o 1.1 1.2 1.3 1.4 2.1 2.2 2.3 - ${outfile}.allElementsBed | sed 's/ /\t/g' | sort -k2b,2 | \
  join -1 2 -2 4 -a1 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 2.1 2.2 2.3 - ${outfile}.allElementsBed | sed 's/ /\t/g' | \
  awk '{dist=int(((\$6+\$7)/2-(\$9+\$10)/2)/1000)}{posdist=dist;if (dist<0) posdist=-dist}{if (posdist!=0) print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"posdist\" kb\"}' > ${outfile}.withGeneNames" >> ${script}
  echo " zcat -f ${outfile}.withGeneNames.tmp | grep Local >> ${outfile}.withGeneNames" >> ${script}

  #and make nicer names for the nodes
  echo "zcat -f ${outfile}.withGeneNames | cut -f1 > ${outfile}.nodes" >> ${script}
  echo "zcat -f ${outfile}.withGeneNames | cut -f2 >> ${outfile}.nodes" >> ${script}
  echo "zcat -f ${outfile}.nodes | sort | uniq | grep -v ":" | awk '{print \$1\"\t\"\$1}' > ${outfile}.nodeLabels" >> ${script}
  echo "zcat -f ${outfile}.nodes | sort | uniq | grep ":" | \
  sed 's/:/\t/g' | sed 's/-/\t/g' | \
  awk '{n=((int(\$2+\$3)/2000)/1000)}{printf \"%s:%s-%s\tchr%s: %.3f kb\n\", \$1,\$2,\$3,\$1,n}' >> ${outfile}.nodeLabels" >> ${script}
  echo "rm ${outfile}.withGeneNames.tmp ${outfile}.TSSsorted"
}

make_detailed_net(){
  local infile=$1
  local outfile=$2
  local pval_thresh=$3
  local LDthresh=$4
  local script=$5

  annos=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits/Annotations/
  peaks2TSS=${annos}/peaks2TSS.txt
  tss=${annos}/pointTSS.bed
  regElements=${annos}/regElts.mergedHistonesDhs.bed

  #testing
  infile=${out}.overlapGWAS.bed
  outfile=${out}.overlapGWAS.network
  pval_thresh=0.00001
  LDthresh=0.8
  script='test.sh'

  #first, GWAS-tagSNP--->QTL-SNP (this is an LD edge)
  echo "zcat -f ${infile} | \
  awk '{if ((\$19<=${pval_thresh}) && (\$9>=${LDthresh})) print \$0}' > ${outfile}.relevant" > ${script}
  echo "zcat -f ${outfile}.relevant | awk '{print \$18\"\t\"\$4\"\tLD\"}' | sort | uniq > ${outfile}.complete" >> ${script}
  #second, QTL----> local reg element
  echo "zcat -f ${outfile}.relevant | awk '{print \$4\"\t\"\$11\"\tLocal\"}' | sort | uniq >> ${outfile}.complete" >> ${script}
  #third, QTL----> distal reg element
  echo "zcat -f ${outfile}.relevant | awk '{print \$4\"\t\"\$12\"\tDistal\"}' | sort | uniq >> ${outfile}.complete" >> ${script}

  #first, GWAS-tagSNP--->QTL-SNP (this is an LD edge)
  echo "zcat -f ${outfile}.relevant | awk '{print \$18\"\t\"\$4\"\tLD\"}' | awk '!seen[$2]++' | sort | uniq > ${outfile}.small" >> ${script}
  #second, QTL----> local reg element
  echo "zcat -f ${outfile}.relevant | awk '{print \$4\"\t\"\$11\"\tLocal\"}' | awk '!seen[$1]++' | sort | uniq >> ${outfile}.small" >> ${script}
  #third, QTL----> distal reg element
  echo "zcat -f ${outfile}.relevant | awk '{print \$4\"\t\"\$12\"\tDistal\"}' | awk '!seen[$1]++' | sort | uniq >> ${outfile}.small" >> ${script}


}

setup_vars ####
columns=chrLDSNP_startLDSNP_endLDSNP_QTL,chr:start-end,LocalDistal,LocalPeak,DistalPeak
GWASPeyton=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-05-30/GWAS_hits_Peyton/
gwases=$(echo $(ls ${GWASPeyton}/expanded_LD_geno/rsq_0.8/*/*bed))
pref=${qtldir}/ALL.LocalandDistalQTLs_expandYRILD
bed_interest=${pref}.2regElts.2TSS.bed
gwas=/srv/gsfs0/projects/snyder/oursu/histoneQTL/GWAS_hits_2015-05-30/GWAS_hits_Peyton//expanded_LD_geno/rsq_0.8/roadmap/EUR.CD_pruned_rsq_0.8_expanded_rsq_0.8.bed


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
 make_net ${out}.overlapGWAS.bed ${out}.overlapGWAS.network 0.00001 0.8 ${s}
 chmod 755 ${s}
 qsub -o ${s}.o -e ${s}.e ${s}
done

#make 1 big GWAS table
for gwas in ${gwases};
do
  out=${gwasdir}/overlapQTL_$(basename ${gwas} | sed 's/.bed//g')/overlapQTL_$(basename ${gwas}).bed
  zcat -f ${out}.overlapGWAS.bed | \
  awk -v gwasname=$(echo $(basename ${gwas}) | sed 's/_no_TAG_SNPs_in_LD_w_SNPs_in_LD.bed//g') '{print $0"\t"gwasname}' > ${out}.GWAStable.txt
done
columns=chrLDSNP_startLDSNP_endLDSNP_QTLSNP,QTLchr:start-end,QTLType,LocalPeak,DistalPeak
echo "${columns},R2,LDSNP,RegEltLocal,RegEltDistal,TSSLocal,TSSDistal_GWASchr_GWASstart_GWASend_GWAStagSNPID_GWASassociationPval_Disease" | \
sed 's/_/\t/g' | sed 's/,/\t/g' > ${gwasdir}/GWAStable.total.txt
zcat -f ${gwasdir}/overlapQTL_*/*GWAStable.txt >> ${gwasdir}/GWAStable.total.txt
zcat -f ${gwasdir}/GWAStable.total.txt | head -n1 > ${gwasdir}/GWAStable.total.GWASsig.txt
zcat -f ${gwasdir}/GWAStable.total.txt | awk '{if ($19<=0.00000001) print $0}' >> ${gwasdir}/GWAStable.total.GWASsig.txt




