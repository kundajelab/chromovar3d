## FUNCTIONS
#############

OverlapRegionsWithTFBS () {
	regionfile=$1
	TFData=$2
	TFData_fileRegex=$3
	outdir=$4
	overlapWindow=$5
	aggregatedOver=$6
	out_prefix=$7
	bashrc=$8
	#make outdir if it doesn't exist
	mkdir -p ${outdir}
	#go through each TF dataset and compute overlap using bedtools
	chips=$(ls ${TFData}/${TFData_fileRegex})
	for chip in ${chips};
   	do
    	chip_base=$(basename ${chip} | sed 's/.gz//g')
    	outname=${outdir}/${out_prefix}.${chip_base}_OVERLAP_$(basename ${aggregatedOver}).tfbsCount

    	s=${outname}.sh
    	echo "source ${bashrc}" > ${s}
    	echo "echo ${chip_base} > ${outname}" >> ${s}
    	#windowbed between our regions and the TF: bedtools window -w ${overlapWindow} -c -a ${regionfile} -b -
    	#keep regions for which we found >=1 TFBS: awk '{if (\$5>0) print \$0}' 
    	#aggregate these over the aggregatedOver bed file: bedtools intersect -c -a ${peakfile} -b - 
    	#keep just the counts: cut -f5 
    	echo "zcat -f ${chip} | sed 's/chr//g' | \
    	bedtools window -w ${overlapWindow} -c -a ${regionfile} -b - | \
    	awk '{if (\$5>0) print \$0}' | \
    	bedtools intersect -c -a ${aggregatedOver} -b - | \
    	cut -f5 >> ${outname}" >> ${s} 
    	chmod 777 ${s}
    	echo $s
    	qsub -o ${s}.o -e ${s}.e ${s}	
    done
}


#Get arguments (just set to NA if not applicable)
regionfile=$1 #where to look for overlaps
TFData=$2 #directory with Tfs for which to compute the overlap
TFData_fileRegex=$3 #TF data name regex
outdir=$4 #directory where to write output
overlapWindow=$5 #window of overlap relative to the $regionfile
aggregatedOver=$6 #aggregate overlaps over this bed file
bashrc=$7 #bashrc that loads bedtools
out_prefix=$8
step=$9 #step: computeOverlaps, overlapsToMatrix, overlapEnrichment
sig_aggregatedOver=${10} #the significant set with which to overlap
bedpe=${11} #for pairwise enrichment

echo "Your step is ${step}"

if [[ $step == "computeOverlaps" ]];then
	echo "======computing overlaps"
	OverlapRegionsWithTFBS ${regionfile} ${TFData} ${TFData_fileRegex} ${outdir} ${overlapWindow} ${aggregatedOver} ${out_prefix} ${bashrc}
fi

overlapMatrix=${outdir}/${out_prefix}.overlapMatrix
if [[ $step == "overlapMatrix" ]];then
	echo "Making the overlap matrix"
	#combine all TFBS overlaps into 1 file
	echo chr_start_end_peak | sed 's/_/\t/g' > ${overlapMatrix}.tmp
	zcat -f ${aggregatedOver} >> ${overlapMatrix}.tmp
	paste ${overlapMatrix}.tmp ${outdir}/${out_prefix}*tfbsCount | gzip > ${overlapMatrix}.gz
	#remove scripts used for the previous step
	rm ${outdir}/${out_prefix}*tfbsCount 
	cat ${outdir}/${out_prefix}*sh > ${outdir}/${out_prefix}.totalScripts 
	rm ${outdir}/${out_prefix}*sh ${outdir}/${out_prefix}*sh.o ${outdir}/${out_prefix}*sh.e ${overlapMatrix}.tmp
fi

#if this is for a pairwise overlap enrichment, then we want to put together the pairwise features starting from a bedpe file
if [[ $step == "mergeOverlapMatrixByBedpe" ]];then 
	s=${outdir}/${out_prefix}.mergeOverlapMatrixByBedpe.sh
	echo "zcat -f ${bedpe} | awk '{print \$1\":\"\$2\"-\"\$3\"\t\"\$0}' | sort -k1b,1 > ${bedpe}.sort1" > ${s}
	#prepare for overlap for interacting site 1
	echo "zcat -f ${outdir}/${out_prefix}.overlapMatrix.gz | head -n1 | awk '{for (i=1;i<=NF;i++) printf \"site1%s \",\$i;printf \"\n\"}' | \
	sed 's/site1chr/chr/g' | sed 's/site1start/start/g' | sed 's/site1end/end/g' > ${outdir}/${out_prefix}.overlapMatrix.gz.site1" >> ${s}
	echo "zcat -f ${outdir}/${out_prefix}.overlapMatrix.gz | grep -v chr >> ${outdir}/${out_prefix}.overlapMatrix.gz.site1" >> ${s}
	echo "zcat -f ${outdir}/${out_prefix}.overlapMatrix.gz.site1 | awk '{print \$1\":\"\$2\"-\"\$3\"\t\"\$0}' | sort -k1b,1 > ${outdir}/${out_prefix}.overlapMatrix.gz.sort1" >> ${s}
	#prepare for overlap for interacting site 2
	echo "zcat -f ${outdir}/${out_prefix}.overlapMatrix.gz | head -n1 | awk '{for (i=1;i<=NF;i++) printf \"site2%s \",\$i;printf \"\n\"}' | \
	sed 's/site2chr/chr/g' | sed 's/site2start/start/g' | sed 's/site2end/end/g' > ${outdir}/${out_prefix}.overlapMatrix.gz.site2" >> ${s}
	echo "zcat -f ${outdir}/${out_prefix}.overlapMatrix.gz | grep -v chr >> ${outdir}/${out_prefix}.overlapMatrix.gz.site2" >> ${s}
	echo "zcat -f ${outdir}/${out_prefix}.overlapMatrix.gz.site2 | awk '{print \$1\":\"\$2\"-\"\$3\"\t\"\$0}' | sort -k1b,1 > ${outdir}/${out_prefix}.overlapMatrix.gz.sort2" >> ${s}
	#join with site 1 data
	echo "join -1 1 -2 1 ${bedpe}.sort1 ${outdir}/${out_prefix}.overlapMatrix.gz.sort1 | sed 's/ /\t/g' | cut --complement -f1 > ${outdir}/${out_prefix}.site1Data" >> ${s}
	echo "zcat -f ${outdir}/${out_prefix}.site1Data | awk '{print \$4\":\"\$5\"-\"\$6\"\t\"\$0}' | sort -k1b,1 > ${outdir}/${out_prefix}.site1Data.sort2" >> ${s}
	#join with site 2 data
	echo "join -1 1 -2 1 ${outdir}/${out_prefix}.site1Data.sort2 ${outdir}/${out_prefix}.overlapMatrix.gz.sort2 | sed 's/ /\t/g' | cut --complement -f1 > ${outdir}/${out_prefix}.site1Data.site2Data" >> ${s}
	#put in the column names first, make a nice file and delete the tmps
	echo "zcat -f ${outdir}/${out_prefix}.site1Data.site2Data | grep chr > ${outdir}/${out_prefix}.site1site2Data" >> ${s}
	echo "zcat -f ${outdir}/${out_prefix}.site1Data.site2Data | grep -v chr >> ${outdir}/${out_prefix}.site1site2Data" >> ${s}
	echo "zcat -f ${outdir}/${out_prefix}.site1site2Data | gzip > ${outdir}/${out_prefix}.site1site2Data.gz" >> ${s}
	echo "rm ${outdir}/${out_prefix}.site1site2Data ${outdir}/${out_prefix}.overlapMatrix.gz.site1 ${outdir}/${out_prefix}.overlapMatrix.gz.site2 \
	${outdir}/${out_prefix}.overlapMatrix.gz.sort1 ${outdir}/${out_prefix}.overlapMatrix.gz.sort2 \
	${outdir}/${out_prefix}.site1Data ${outdir}/${out_prefix}.site1Data.sort2 ${outdir}/${out_prefix}.site1Data.site2Data" >> ${s}
	chmod 755 ${s}
	qsub -o ${s}.o -e ${s}.e ${s}
fi

if [[ $step == 'enrichmentAnalysisByBedpe' ]];then
	tfs=$(zcat -f ${outdir}/${out_prefix}.site1site2Data.gz  | head -n1 | sed 's/\t/\n/g' | grep site1 | grep -v peak | sed 's/site1//g')
	echo $(zcat -f ${outdir}/${out_prefix}.site1site2Data.gz  | head -n1 | sed 's/\t/\n/g' | grep site1 | grep -v peak | wc -l)
	for tf1 in ${tfs};
	do 
		s=${outdir}/${out_prefix}.enrichmentAnalysisByBedpe.${tf1}.sh
		for tf2 in ${tfs};
		do
			
		 	#tf1=MergedPeaks_ChIPseq_CTCF
			#tf2=MergedPeaks_ChIPseq_SPI1
			tf1_1=$(zcat -f ${outdir}/${out_prefix}.site1site2Data.gz  | head -n1 | sed 's/\t/\n/g' | grep -n site1${tf1} | cut -f1 -d ":")
			tf1_2=$(zcat -f ${outdir}/${out_prefix}.site1site2Data.gz  | head -n1 | sed 's/\t/\n/g' | grep -n site2${tf1} | cut -f1 -d ":")
			tf2_1=$(zcat -f ${outdir}/${out_prefix}.site1site2Data.gz  | head -n1 | sed 's/\t/\n/g' | grep -n site1${tf2} | cut -f1 -d ":")
			tf2_2=$(zcat -f ${outdir}/${out_prefix}.site1site2Data.gz  | head -n1 | sed 's/\t/\n/g' | grep -n site2${tf2} | cut -f1 -d ":")
			echo "zcat -f ${outdir}/${out_prefix}.site1site2Data.gz | \
			awk -F \"\t\" '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6\"\t\"\$${tf1_1}\"\t\"\$${tf1_2}\"\t\"\$${tf2_1}\"\t\"\$${tf2_2}}' | \
			awk '{if ((\$7>0 && \$10>0) || (\$8>0 && \$9>0)) print \$0\"\tpair\";else print \$0\"\tNotPair\"}' | gzip > ${outdir}/${out_prefix}.${tf1}.${tf2}.gz" >> ${s}
			# make a file with the following quantities
			# total interactions, freq(TF1_1), freq (TF1_2), freq (TF2_1), freq (TF2_2), pair(1,2), pair(2,1)
			echo "f11=\$(zcat -f ${outdir}/${out_prefix}.${tf1}.${tf2}.gz | awk '{if (\$7>0) print \$7}' | wc -l)" >> ${s}
			echo "f12=\$(zcat -f ${outdir}/${out_prefix}.${tf1}.${tf2}.gz | awk '{if (\$8>0) print \$8}' | wc -l)" >> ${s}
			echo "f21=\$(zcat -f ${outdir}/${out_prefix}.${tf1}.${tf2}.gz | awk '{if (\$9>0) print \$9}' | wc -l)" >> ${s}
			echo "f22=\$(zcat -f ${outdir}/${out_prefix}.${tf1}.${tf2}.gz | awk '{if (\$10>0) print \$10}' | wc -l)" >> ${s}
			echo "p12=\$(zcat -f ${outdir}/${out_prefix}.${tf1}.${tf2}.gz | awk '{if ((\$7>0 && \$10>0)) print \$0}' | wc -l)" >> ${s}
			echo "p21=\$(zcat -f ${outdir}/${out_prefix}.${tf1}.${tf2}.gz | awk '{if ((\$8>0 && \$9>0)) print \$0}' | wc -l)" >> ${s}
			echo "total=\$(zcat -f ${outdir}/${out_prefix}.${tf1}.${tf2}.gz | wc -l)" >> ${s}
			echo "echo Item1_Item2_Total_F1atsite1_F1atsite2_F2atsite1_F2atsite2_Paired1atsite1with2atsite2_Paired1atsite2with2atsite1 | sed 's/_/\t/g' > ${outdir}/${out_prefix}.${tf1}.${tf2}.summary.txt" >> ${s}
			echo "echo ${tf1}__${tf2}__\${total}__\${f11}__\${f12}__\${f21}__\${f22}__\${p12}__\${p21} | sed 's/__/\t/g' >> ${outdir}/${out_prefix}.${tf1}.${tf2}.summary.txt" >> ${s}
		done
		chmod 755 ${s}
		qsub -o ${s}.o -e ${s}.e ${s}
	done
fi

if [[ $step == "enrichmentAnalysis" ]];then
	source ${bashrc}
	#perform enrichment analysis
	echo "computing significance"
	${R_WITH_PACKAGES} ${enrichRcode} ${overlapMatrix}.gz ${sig_aggregatedOver} ${outdir}/${out_prefix}.overlapEnrichIN$(basename ${sig_aggregatedOver})
fi


