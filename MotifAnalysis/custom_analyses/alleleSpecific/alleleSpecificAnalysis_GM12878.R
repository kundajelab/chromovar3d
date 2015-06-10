


infile=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Local/SNPQTLmatrix/
TFdata=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/TFBS/alleleSpecificChIP/NA12878_AS_SNPs.hg19.vcf.gz_allTFs
cat ${TFdata} | sort -k3b,3 > ${TFdata}.sort
for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA
do
 zcat ${infile}/SNPQTLmatrix.${hmark}.gz | sort -k3b,3 > ${infile}/SNPQTLmatrix.${hmark}.gz.sort
 zcat ${infile}/SNPQTLmatrix.${hmark}.gz | head -n1 | awk '{print $0"\tchr\tSNPpos\tSNPID\tTF\tREF\tALT"}' > ${TFdata}.with${hmark}SNPs
 join -1 3 -2 3 ${infile}/SNPQTLmatrix.${hmark}.gz.sort ${TFdata}.sort | sed 's/ /\t/g' >> ${TFdata}.with${hmark}SNPs
 rm ${TFdata}.with${hmark}SNPs.gz
 gzip ${TFdata}.with${hmark}SNPs
done

#count peaks explained
for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA
do
 echo ${hmark}
 echo "Number of peaks explained for ${hmark} $(zcat ${TFdata}.with${hmark}SNPs.gz | cut -f4 | sort | wc -l)" > ${TFdata}.with${hmark}SNPs.peakCount
 samedir=$(zcat ${TFdata}.with${hmark}SNPs.gz | awk '{samedir=$5*($16-$15)}{print $0"\t"samedir}' | grep -v ALT | awk '{if ($17>=0) print $17}' | wc -l)
 diffdir=$(zcat ${TFdata}.with${hmark}SNPs.gz | awk '{samedir=$5*($16-$15)}{print $0"\t"samedir}' | grep -v ALT | awk '{if ($17<0) print $17}' | wc -l)
 echo "Same direction for ${hmark}: ${samedir}" > ${TFdata}.with${hmark}SNPs.directionCount
 echo "Diff direction for ${hmark}: ${diffdir}" >> ${TFdata}.with${hmark}SNPs.directionCount
 echo "Percent same direction: $(echo "nothing" | awk -v s=${samedir} -v d=${diffdir} '{print s/(s+d)}')" >> ${TFdata}.with${hmark}SNPs.directionCount
done

#cat ${TFdata}.with*SNPs.peakCount
#Number of peaks explained for dhs 1502
#Number of peaks explained for H3K27AC 565
#Number of peaks explained for H3K4ME1 887
#Number of peaks explained for H3K4ME3 385
#Number of peaks explained for RNA 106

#direction of TF binding (counting snp-TF pairs): 
#cat ${TFdata}.with*SNPs.directionCount
#Same direction for dhs: 862
#Diff direction for dhs: 639
#Percent same direction: 0.574284

#Same direction for H3K27AC: 313
#Diff direction for H3K27AC: 251
#Percent same direction: 0.554965

#Same direction for H3K4ME1: 519
#Diff direction for H3K4ME1: 367
#Percent same direction: 0.585779

#Same direction for H3K4ME3: 205
#Diff direction for H3K4ME3: 179
#Percent same direction: 0.533854

#Same direction for RNA: 51
#Diff direction for RNA: 54
#Percent same direction: 0.485714
