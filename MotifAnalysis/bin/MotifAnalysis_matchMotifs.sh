
#!/bin/bash
# ========================================
# Read in arguments and check for errors
# ========================================
if [[ "$#" -lt 6 ]]
    then
    echo "Usage: $(basename $0)">&2
    echo '<PWM> PWM name'>&2
    echo '<PWMtemplate> Template of pwm path'>&2
    echo '<source place>'>&2
    echo '<out>'>&2
    echo '<regions>'>&2
    echo '<chips>'>&2
    exit 1
fi

PWM=$1
PWM_template=$2
source_place=$3
OUT=$4
source ${source_place}
regions=$5
chips=$6

#Output files
mkdir -p ${OUT}

echo "Starting"
echo ${OUT}/MotifMatch_${PWM}.bed.gz

#Get genomic matches to the motif
homermatch=${OUT}/MotifMatch_${PWM}.bed.gz
if [ ! -e "${homermatch}" ]; then
    matchScript=${OUT}/MotifMatch_${PWM}.bed.gz.sh
    echo ${matchScript}
    echo "source ${source_place}" > ${matchScript}
    echo "export PATH="'${HOMERDIR}'":"'$PATH' >> ${matchScript}
    tempo=/srv/gsfs0/scratch/oursu/
    echo "export TMP=${tempo}" >> ${matchScript}
    echo "cd "'${TMP}' >> ${matchScript}
    #use the PWM with match threshold=0                                                                                                                                     
    PWM_for_match=$(echo ${PWM_template} | sed 's/PWMNAME/'${PWM}'/g')
    echo "zcat ${regions} | "'${HOMERDIR}'"/annotatePeaks.pl - hg19 -noann -nogene -size given -m ${PWM_for_match} -mbed ${homermatch}_notgz" >> ${matchScript}
    echo "cat ${homermatch}_notgz | gzip > ${homermatch}" >> ${matchScript}
    echo "zcat ${homermatch} | bedtools intersect -u -a - -b ${chips} | gzip > ${homermatch}.OverlapChIPseq.gz" >> ${matchScript}
    echo "rm ${matchScript}.o ${homermatch}_notgz ${matchScript}.e" >> ${matchScript}
    echo "rm ${chips}" >> ${matchScript}
    qsub -l h_vmem=10G -o ${matchScript}.o -e ${matchScript}.e ${matchScript} 
fi


