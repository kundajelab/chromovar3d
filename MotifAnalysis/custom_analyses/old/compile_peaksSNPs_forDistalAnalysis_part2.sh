
distaldata=/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-07-28/HiC_local_QTL_dist_QTL.all.txt_augmentPeak1Peak2IDs.ChrGenePeak2.PlusOverlappedRNA.txt
cat ${distaldata} | sort -k3 > ${distaldata}.sorted
for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs
do
    zcat /srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-07-28/SNPQTLmatrix.${hmark}.gz | awk -v histmark=$hmark '{print histmark"_"$4"\t"$2"\t"$3"\t"$1}' | sort -k1b,1 | gzip > /srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-07-28/SNPQTLmatrix.${hmark}.sorted.gz
    echo "peak2id_chr_SNP_snp.position_chromo_gene" | tr "_" "\t" > /srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-07-28/distal/SNPQTLmatrix.${hmark}
    zcat /srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-07-28/SNPQTLmatrix.${hmark}.sorted.gz | join -1 1 -2 3 - ${distaldata}.sorted >> /srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-07-28/distal/SNPQTLmatrix.${hmark}
    cat /srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-07-28/distal/SNPQTLmatrix.${hmark} | gzip > /srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-07-28/distal/SNPQTLmatrix.${hmark}.gz
    rm /srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-07-28/distal/SNPQTLmatrix.${hmark}
    rm /srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/results/2014-07-28/SNPQTLmatrix.${hmark}.sorted.gz
done
rm ${distaldata}.sorted
