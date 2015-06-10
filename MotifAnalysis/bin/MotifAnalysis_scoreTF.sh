
#!/bin/bash
# ========================================
# Read in arguments and check for errors
# ========================================
if [[ "$#" -lt 8 ]]
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

#Output files
mkdir ${OUT}
outsuf=$(basename ${OUT})
motifScoreDir=${OUT}/$(basename ${PWM})
mkdir ${motifScoreDir}
mkdir ${motifScoreDir}/Scores
scriptName=${motifScoreDir}/Scores/$(basename ${PWM})_scoreTF

#Figure out motif length
PWMLen=`cat ${PWM} | wc -l`
motifLen=`expr $PWMLen - 1`

# ----- Score the motif
#----------------------
#Remember the regions file
out=${motifData}/SNP2bed_motifLen${motifLen}/$(basename ${SNPQTL})_motifLen${motifLen}
past_output=${out}.snpsite
for filter in $(echo ${filters} | sed -n 1'p' | tr ',' '\n')
do
    echo ${filter}
    cur_output=${past_output}.filter.$(echo $(basename ${filter}) | sed 's/.bed.gz//g')
    past_output=${cur_output}
done
regions_file=${cur_output}.bed.gz

#Score each individual's 2 genomes	    
#=================================
compiled_maternal=${motifScoreDir}/Scores/$(basename ${PWM})_maternalScores.compiled
compiled_paternal=${motifScoreDir}/Scores/$(basename ${PWM})_paternalScores.compiled
echo "Scoring TF ${PWM}"
counter=1
tempo=${motifData}_tempo
mkdir ${tempo}
while read line
do      

    items=($(echo ${line} | tr ' ' '\n'))
    person=${items[0]}
    scoring_cur=${scriptName}_${person}.sh
    rm ${scoring_cur}.e
    rm ${scoring_cur}.o
    echo "source ${source_place}" > ${scoring_cur}
    person_fasta=${motifData}/SNP2bed_motifLen${motifLen}/Fasta/${person}
    basename_regions=$(basename ${regions_file})
    maternal_fa=${person_fasta}.maternalSeqsIN.${basename_regions}.fa
    paternal_fa=${person_fasta}.paternalSeqsIN.${basename_regions}.fa
    maternal_out=${motifScoreDir}/Scores/${person}_maternalScores_$(basename ${PWM}).gz
    paternal_out=${motifScoreDir}/Scores/${person}_paternalScores_$(basename ${PWM}).gz
    ##Here start the commands
    echo "export PATH="'${HOMERDIR}'":"'$PATH' >> ${scoring_cur}
    #for regions in ${regions_tf}part*;
    #do
    #tempo=${compiled_maternal}_${person}
    #mkdir ${tempo}
    tempo=/srv/gsfs0/scratch/oursu/${person}_$(basename ${PWM})
    mkdir ${tempo}
    echo "export TMP=${tempo}"  >> ${scoring_cur}
    echo "cd "'${TMP}' >> ${scoring_cur}
    echo '${HOMERDIR}'"/homer2 find -i ${maternal_fa} -mscore -m ${PWM} | gzip > ${maternal_out}" >> ${scoring_cur}
    echo '${HOMERDIR}'"/homer2 find -i ${paternal_fa} -mscore -m ${PWM} | gzip > ${paternal_out}" >> ${scoring_cur}
    echo "lines_maternal="'`zcat '"${maternal_out}"' | wc -l`' >> ${scoring_cur}
    echo "if [[ "'$lines_maternal'" -eq 0 ]]; then "'${HOMERDIR}'"/homer2 find -i ${maternal_fa} -mscore -m ${PWM} | gzip > ${maternal_out};fi" >> ${scoring_cur} 
    echo "lines_paternal="'`zcat '"${paternal_out}"' | wc -l`' >> ${scoring_cur}
    echo "if [[ "'$lines_paternal'" -eq 0 ]]; then "'${HOMERDIR}'"/homer2 find -i ${paternal_fa} -mscore -m ${PWM} | gzip > ${paternal_out};fi" >> ${scoring_cur}
    echo "rm -r ${tempo}" >> ${scoring_cur}
    #echo " cat ${scoring_cur}.e ${scoring}.o > ${scoring_cur}.o" >> ${scoring_cur}
    #echo " cat ${scoring_cur}.e > ${scoring_cur}.o" >> ${scoring_cur}
    #echo "rm ${scoring_cur}.e ${scoring_cur}.o ${scoring_cur}" >> ${scoring_cur}
    #done
    #echo "cat ${maternal_out}.notgz | gzip > ${maternal_out}" >> ${scoring_cur}
    #echo "cat ${paternal_out}.notgz | gzip > ${paternal_out}" >> ${scoring_cur}
    counter=`expr ${counter} + 1`
    chmod 711 ${scoring_cur}
    #qsub -l hostname=scg3-0-* -o ${scoring_cur}.o -e ${scoring_cur}.e ${scoring_cur}
    qsub -wd /srv/gsfs0/scratch/oursu/ -l hostname=scg3-0-* -o ${scoring_cur}.o -e ${scoring_cur}.e ${scoring_cur}
done < ${metadata}




