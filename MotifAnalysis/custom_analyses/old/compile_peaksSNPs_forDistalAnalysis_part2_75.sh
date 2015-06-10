
distaldata=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-12/DistalQTL.txt.ChrGenePeak2.PlusOverlappedRNA.txt
cat ${distaldata} | sort -k3 > ${distaldata}.sorted
DATA=/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-12/
mkdir ${DATA}distal
for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA
do
    zcat ${DATA}SNPQTLmatrix.${hmark}.gz | awk -v histmark=$hmark '{print histmark"_"$4"\t"$2"\t"$3"\t"$1}' | sort -k1b,1 | gzip > ${DATA}/SNPQTLmatrix.${hmark}.sorted.gz
    echo "gene_chr_SNP_snp.position_chromo_localpeak" | tr "_" "\t" > ${DATA}/distal/SNPQTLmatrix.${hmark}
    zcat ${DATA}/SNPQTLmatrix.${hmark}.sorted.gz | join -1 1 -2 3 - ${distaldata}.sorted >> ${DATA}/distal/SNPQTLmatrix.${hmark}
    cat ${DATA}/distal/SNPQTLmatrix.${hmark} | gzip > ${DATA}/distal/SNPQTLmatrix.${hmark}.gz
    rm ${DATA}/distal/SNPQTLmatrix.${hmark}
    rm ${DATA}/SNPQTLmatrix.${hmark}.sorted.gz
done
rm ${distaldata}.sorted
