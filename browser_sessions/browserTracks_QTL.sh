module load tabix

#==== Peak files ======
#======================
outmain=/srv/gs1/projects/snyder/jzaugg/histoneQTL/chromovar3d
#Prepare peak files
for hmark in H3K4ME1 H3K4ME3 H3K27AC 
do
 out=${outmain}/histone/results/peaks/merged_peaks/${hmark}_mergedPeaks.bed
 zcat -f ${outmain}/histone/results/peaks/merged_peaks/mergePeaks_${hmark}.gzremoveBlacklist.gz_Log10PvalueThreshold_5.gz | \
 awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"NR"\t."}' > ${out}.tmp
 cat ${out}.tmp | sort -k1,1 -k2n,2n -k3n,3n | awk '{print $0"\t"NR"\t."}' > ${out}
 bgzip -c ${out} > ${out}.gz
 tabix -p bed ${out}.gz
done
#for dhs
out=${outmain}/DNase/data/hg19/peaks/subsampled/mergedPeaks_TrimFromSummit/DNase_mergedPeaks.bed
zcat -f ${outmain}/DNase/data/hg19/peaks/subsampled/mergedPeaks_TrimFromSummit/DNase_TrimFromSummit.mergeBedremoveBlacklist_Log10PvalueThreshold_5.gz | \
awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"NR"\t."}' > ${out}.tmp
cat ${out}.tmp | sort -k1,1 -k2n,2n -k3n,3n | awk '{print $0"\t"NR"\t."}' > ${out}
bgzip -c ${out} > ${out}.gz
tabix -p bed ${out}.gz

#==== Pairwise files =========
#=============================
pairwise_track(){
  local snpqtl_local=$1
  local peak_anno=$2
  local out=$3
  zcat -f ${peak_anno} | sort -k4b,4 > ${peak_anno}_unzip
  zcat ${snpqtl_local} | awk -F "\t" '{if ($9=="pass") print $10"_"$4"\t"$0}' | \
   sort -k1b,1 | join -a1 -1 1 -2 4 - ${peak_anno}_unzip | sed 's/ /\t/g' | \
    awk -F "\t" '{minusone=$2-1}{print "chr"$3"\t"minusone"\t"$2"\tchr"$3":"$14"-"$15","$6}' > ${out}.tmp
  zcat -f ${peak_anno} | sort -k4b,4 > ${peak_anno}_unzip
  zcat ${snpqtl_local} | awk -F "\t" '{if ($9=="pass") print $10"_"$4"\t"$0}' | \
   sort -k1b,1 | join -a1 -1 1 -2 4 - ${peak_anno}_unzip | sed 's/ /\t/g' | \
    awk -F "\t" '{minusone=$2-1}{print "chr"$3"\t"$14"\t"$15"\tchr"$3":"minusone"-"$2","$6}' >> ${out}.tmp
    cat ${out}.tmp | sort -k1,1 -k2n,2n -k3n,3n | awk '{print $0"\t"NR"\t."}' > ${out}
    bgzip -c ${out} > ${out}.gz
    tabix -p bed ${out}.gz
}

#Local QTLs
#===========
outmain=/srv/gs1/projects/snyder/jzaugg/histoneQTL/chromovar3d
for hmark in H3K4ME1 H3K4ME3 H3K27AC dhs RNA
do
  echo ${hmark}
  out=/srv/gsfs0/projects/snyder/oursu/histoneQTL/browserTracks/2015-07-17/pairwise/QTL/LocalQTLs.${hmark}
  pairwise_track ${outmain}/motifs/results/2015-06-13/Local/SNPQTLmatrix/SNPQTLmatrix.${hmark}.gz ${outmain}/motifs/results/2015-06-13/OverlapEnrichment/Peaks/peaks.${hmark}.anno.gz ${out}
done

#Distal QTLs
#===========
/home/oursu/devtools/R-3.0.2/bin/Rscript /srv/gsfs0/projects/kundaje/users/oursu/code/ChromatinVariationCode/browser/distalQTL_to_pairwiseTrack.R
distals=/srv/gsfs0/projects/snyder/oursu/histoneQTL/browserTracks/2015-07-17/pairwise/QTL/
distpre=DistalQTL_2015-05-30.states.prox.frag.FDR10.allChr.2000kb.HiC_TH0.4.LDcorrected_0.2.states.BROWSER.
for hmark in H3K4ME3 #H3K4ME1 H3K4ME3 H3K27AC dhs RNA
do
 out=${distals}/DistalQTLs.${hmark}
 zcat -f ${distals}/${distpre}${hmark}.pairwiseTrack | sed 's/:/\t/g' | sed 's/-/\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$4":"$5"-"$6","$7}' > ${out}.tmp
 zcat -f ${distals}/${distpre}${hmark}.pairwiseTrack | sed 's/:/\t/g' | sed 's/-/\t/g'| awk '{print $4"\t"$5"\t"$6"\t"$1":"$2"-"$3","$7}' >> ${out}.tmp
 cat ${out}.tmp | sort -k1,1 -k2n,2n -k3n,3n | awk '{print $0"\t"NR"\t."}' > ${out}
 bgzip -c ${out} > ${out}.gz
 tabix -p bed ${out}.gz
done

#copy the above to the browserFiles directory

#ChIA-PET
#========
bedpe_to_track (){
 local bedpe=$1
 out=${bedpe}.pairwiseTrack

 zcat -f ${bedpe} | awk '{print $1"\t"$2"\t"$3"\t"$4":"$5"-"$6",1"}' > ${out}.tmp
 zcat -f ${bedpe} | awk '{print $4"\t"$5"\t"$6"\t"$1":"$2"-"$3",1"}' >> ${out}.tmp

 module load tabix
 cat ${out}.tmp | sort -k1,1 -k2n,2n -k3n,3n | awk '{print $0"\t"NR"\t."}' > ${out}
 bgzip -c ${out} > ${out}.gz
 tabix -f -p bed ${out}.gz
 rm ${out}.tmp
}

c1=/srv/gs1/projects/snyder/jzaugg/histoneQTL/chromovar3d/ChIA-PET/mango/FG.cell_H3K4Me3_GM78_std_0.0.interactions.fdr-2.bedpe
c2=/srv/gs1/projects/snyder/jzaugg/histoneQTL/chromovar3d/ChIA-PET/mango/FG.cell_RAD21_GM78_std_0.0.interactions.fdr-2.bedpe
bedpe_to_track ${c1}
bedpe_to_track ${c2}
cp /srv/gs1/projects/snyder/jzaugg/histoneQTL/chromovar3d/ChIA-PET/mango/*gz* /srv/gs1/projects/snyder/jzaugg/histoneQTL/chromovar3d/browserTracks/ChIAPET/