
#!/bin/bash
# ========================================
# Read in arguments and check for errors
# ========================================
if [[ "$#" -lt 10 ]]
    then
    echo "Usage: $(basename $0)">&2
    echo '<PWM> File with motifPWM.'>&2
    echo '<genomic regions for filtering> '>&2
    echo '<metadata file> First column specifying which individuals to use in the analysis. 1. Individual. 2. Maternal genome. 3. Paternal genome.'>&2
    echo '<SNP-QTL> file, gzipped'>&2
    echo '<motifData> directory where to write SNP length files and fastas for the people'>&2
    echo '<source place>'>&2
    echo '<out>'>&2
    echo '<signalfile>' >&2
    echo '<motifMatchFilterBed> File with motif matches across the genome, e.g. TFBS' >&2
    echo '<distalLocal> Prefix of correlation files' >&2
    exit 1
fi

PWM=$1
filters=$2
metadata=$3
SNPQTL=$4
motifData=$5
source_place=$6
OUT=$7
signalfile=$8
motifMatchFilterBed=$9
distallocal=${10}
SNPQTLmatched=${11}
motifMAF=0.05

#Output files
mkdir ${OUT}
outsuf=$(basename ${OUT})
motifScoreDir=${OUT}/$(basename ${PWM})
mkdir ${motifScoreDir}
mkdir ${motifScoreDir}/Scores
scriptName=${motifScoreDir}/Scores/$(basename ${PWM})
mkdir ${motifScoreDir}/Correlations

#Figure out motif length
PWMLen=`cat ${PWM} | wc -l`
motifLen=`expr $PWMLen - 1`

#Remember the regions file
#out=${motifData}/SNP2bed_motifLen${motifLen}/$(basename ${SNPQTL})_motifLen${motifLen}
#past_output=${out}.snpsite
#for filter in $(echo ${filters} | sed -n 1'p' | tr ',' '\n')
#do
#    echo ${filter}
#    cur_output=${past_output}.filter.$(echo $(basename ${filter}) | sed 's/.bed.gz//g')
#    past_output=${cur_output}
#done
#regions_file=${cur_output}.bed

#### Compile motif scores only if they don't exist already compiled
#basename_regions=$(basename ${regions_file})
echo "compiled regions"
echo ${motifScoreDir}/Scores/$(basename ${PWM})_maternalScores.compiled
compiled_maternal=$(echo ${motifScoreDir}/Scores/$(basename ${PWM})_maternalScores.compiled | sed 's/LocalPeakIs//g')
compiled_paternal=$(echo ${motifScoreDir}/Scores/$(basename ${PWM})_paternalScores.compiled | sed 's/LocalPeakIs//g')
if [ ! -e "${compiled_maternal}.gz" ]; then
    #Merging part
    firstTime=true
    scoring_script=${scriptName}_mergeScores.sh
    echo "source ${source_place}" > ${scoring_script}
    while read line
    do
	echo "reading"
	items=($(echo ${line} | tr ' ' '\n'))
	person=${items[0]}
	echo "rm ${scriptName}_motifScoring_${person}.sh" >> ${scoring_script}
	maternal_out=${motifScoreDir}/Scores/${person}_maternalScores_$(basename ${PWM}).gz
	paternal_out=${motifScoreDir}/Scores/${person}_paternalScores_$(basename ${PWM}).gz
	cur_maternal=${maternal_out}.cur
	cur_paternal=${paternal_out}.cur
	echo "echo \"aSNPregionDELIMITE${person}\" | sed 's/DELIMITE/\t/g' > ${cur_maternal}">> ${scoring_script}
	echo "zcat ${maternal_out} | awk '{print "'$1'"\"\t\""'$NF'"}' | sort >> ${cur_maternal}" >> ${scoring_script}
	echo "echo \"aSNPregionDELIMITE${person}\" | sed 's/DELIMITE/\t/g' > ${cur_paternal}">> ${scoring_script}
	echo "zcat ${paternal_out} | awk '{print "'$1'"\"\t\""'$NF'"}' | sort >> ${cur_paternal}" >> ${scoring_script}
	if [ "$firstTime" = true ] ; then
            echo "cat ${cur_maternal} | awk '{print "'$1'"}' | sort | uniq > ${compiled_maternal}" >> ${scoring_script}
            echo "cat ${cur_paternal} | awk '{print "'$1'"}' | sort | uniq > ${compiled_paternal}" >> ${scoring_script}
            firstTime=false
	fi
	echo "join -1 1 -2 1 ${compiled_maternal} ${cur_maternal} > ${compiled_maternal}new" >> ${scoring_script}
	echo "mv ${compiled_maternal}new ${compiled_maternal}" >> ${scoring_script}
	echo "rm ${cur_maternal}" >> ${scoring_script}
	echo "join -1 1 -2 1 ${compiled_paternal} ${cur_paternal} > ${compiled_paternal}new" >> ${scoring_script}
	echo "mv ${compiled_paternal}new ${compiled_paternal}" >> ${scoring_script}
	echo "rm ${cur_paternal}" >> ${scoring_script}
    done < ${metadata}
    echo "cat ${compiled_maternal} | gzip > ${compiled_maternal}.gz" >> ${scoring_script}
    echo "cat ${compiled_paternal} | gzip > ${compiled_paternal}.gz" >> ${scoring_script}
    echo "rm ${compiled_maternal} ${compiled_paternal}" >> ${scoring_script}
    #echo "rm ${motifScoreDir}/Scores/*_maternalScores_$(basename ${PWM}).gz" >> ${scoring_script}
    #echo "rm ${motifScoreDir}/Scores/*_paternalScores_$(basename ${PWM}).gz" >> ${scoring_script}
    #echo "rm ${motifScoreDir}/Scores/*scoreTF_NA*.sh.o" >> ${scoring_script}
    rm ${scoring_script}.e
    rm ${scoring_script}.o
    #Submit first part here                                                                                                                                                      
    merge_jobID=`qsub -l hostname=scg3-0-* -o ${scoring_script}.o -e ${scoring_script}.e ${scoring_script} | cut -d ' ' -f 3 `
    echo "Compiling TF scores because compiled scores do not already exist"
    waitmsg='-hold_jid '${merge_jobID}
fi

#Run significance analysis. Compute true correlations and permutations.
results=${motifScoreDir}/Correlations
# ------ Significance                                                                                                                                                        
outtrue=${results}/$(basename ${PWM}).motifMAF${motifMAF}.correlation.${distallocal}.tab
sig_script=${outtrue}_computeAssociations_${distallocal}.sh
echo "source ${source_place}" > ${sig_script}
echo '${R_WITH_PACKAGES}'" "'${CODEDIR}'"r_code/computeAssociations.R ${metadata} ${PWM} ${SNPQTL} ${compiled_maternal}.gz,${compiled_paternal}.gz ${signalfile} ${outtrue} ${motifMatchFilterBed} ${motifMAF}" >> ${sig_script}
echo "cat ${outtrue}permCorrTable.txt | gzip > ${outtrue}permCorrTable.gz" >> ${sig_script}
echo "rm ${outtrue}permCorrTable.txt"  >> ${sig_script}
qsub ${waitmsg} -l h_vmem=20G -o ${sig_script}.o -e ${sig_script}.e ${sig_script}

