#1. Download data
/srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/src/download_mapped_read_call.sh 
#Output: /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/Degner_mapped_reads/eqtl.uchicago.edu/dsQTL_data/MAPPED_READS/

#2. Lift over from hg18 to hg19
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/src/processDNase.py --step_to_perform liftOver
#Output: /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/liftOver/

#3. Shift reads to use for MACS
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/src/processDNase.py --step_to_perform shiftReads
#Output:/srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/shiftedReads/

#4. Subsample to 35M reads
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/src/processDNase.py --step_to_perform subsample --subsample_reads 35
#Output: /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/subsampled/

#5. Call peaks using MACS2
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/src/processDNase.py --step_to_perform MACS --peaks_on_subsampled --subsample_reads 35
#Output: /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/peaks/subsampled/

#6. Get raw signal
source /srv/gsfs0/projects/kundaje/commonRepository/src/lab_bashrc
activateMCR
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/src/processDNase.py --peaks_on_subsampled --step_to_perform align2rawsignal --subsample_reads 35
#Output: /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/subsampled/align2rawsignal/

#7. Merge peaks
#Trim peaks
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/src/processDNase.py --step_to_perform TrimFromSummit --peaks_on_subsampled --subsample_reads 35
#Output: /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/peaks/subsampled/
#merge
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/src/processDNase.py --step_to_perform mergePeaks --peaks_on_subsampled --subsample_reads 35 --merge_type TrimFromSummit
#Output: /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/peaks/subsampled/mergedPeaks_TrimFromSummit/DNase_TrimFromSummit.mergeBedremoveBlacklist_Log10PvalueThreshold_5.gz

#8. Extract signal
source /srv/gsfs0/projects/kundaje/commonRepository/src/lab_bashrc
activateMCR
python /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/src/processDNase.py --step_to_perform extractSignal --peaks_on_subsampled --subsample_reads 35 --mergedPeaks /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/peaks/subsampled/mergedPeaks_TrimFromSummit/DNase_TrimFromSummit.mergeBedremoveBlacklist_Log10PvalueThreshold_5.gz
#Output: /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/subsampled/align2rawsignal/extractSignal/

#9. Make big data matrix with extracted signal inside merged peaks
#Output:/srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/subsampled//align2rawsignal/extractSignal/DNase_removeBlacklist_Log10PvalueThreshold_5_DATA_MATRIX

datastart=/srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/subsampled//align2rawsignal/extractSignal
cd ${datastart}
dataend=.subsampleTo35000000.wiggler.norm5.rawsignal.mat_VS_mergePeaks_DNase_TrimFromSummit.mergeBedremoveBlacklist_Log10PvalueThreshold_5.MEAN.cagt
for file in *${dataend}; do echo "nameq${file}" | tr 'q' '\t' >tempfile; cat $file | sort | awk '{print $1"_"$2"_"$3"\t"$4}' >>tempfile; mv tempfile ${file}_signal; done

peaks=/srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/data/hg19/peaks/subsampled/mergedPeaks_TrimFromSummit/DNase_TrimFromSummit.mergeBedremoveBlacklist_Log10PvalueThreshold_5.gz
echo "name_chr_start_end" | tr '_' '\t' > ${datastart}/DNaseremoveBlacklist_peaks
zcat ${peaks} | awk '{a=$2+1}{print $1"_"a"_"$3"\t"$1"\t"a"\t"$3}' | sort >> ${datastart}/DNaseremoveBlacklist_peaks

final=${datastart}/DNase_removeBlacklist_Log10PvalueThreshold_5_DATA_MATRIX
initial=${final}_long
cat ${datastart}/DNaseremoveBlacklist_peaks > ${initial}
for file in *${dataend}_signal; do join ${initial} ${file} > ${initial}_new; mv ${initial}_new ${initial};done 
cat ${initial} | perl -pe 's/Rep_0_//g' | perl -pe 's/.subsampleTo35000000.wiggler.norm5.rawsignal.mat_VS_mergePeaks_DNase_TrimFromSummit.mergeBedremoveBlacklist_Log10PvalueThreshold_5.MEAN.cagt//g' | cut -d " " -f2- > ${final}
rm ${final}_long