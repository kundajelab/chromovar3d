#!/bin/bash
# ========================================
# Read in arguments and check for errors
# ========================================
if [[ "$#" -lt 6 ]]
    then
    echo "Usage: $(basename $0)">&2
    echo '<motifLens> Motif lengths separated by comma. This is a list of all lengths of motifs you expect in your file. Example: 8,9,10,11,12'>&2
    echo '<genomic regions for filtering> Bed file. This is a genomic region with which to filter your SNPs. For instance, if your SNPs cover the whole genome, this filter may represent the set of active enhancers in your cell type of interest. This code wil only keep the SNPs that fall inside these filter genomic regions. You can specify multiple bed files, separated by comma, and you wil get the intersection between them.'>&2
    echo '<metadata file> First column specifying which individuals to use in the analysis. 1. Individual. 2. Maternal genome. 3. Paternal genome. You can specify the reference genome for both, if you do not have personal genomes. The paths to a genome is the path to ONR fasta file containing all chromosomes of the genome. '>&2
    echo '<SNP-QTL> file, gzipped. Look at my example.'>&2
    echo '<motifData> directory where to write SNP length files and fastas for the people'>&2
    echo '<source place> bashrc file'>&2
    exit 1
fi

motifLens=$1
filters=$2
metadata=$3
SNPQTL=$4
motifData=$5
source_place=$6

#----- SNP2bed. Now, take all SNPs in the provided file, and make a bed file with the regions of motifLen around the SNP
mkdir ${motifData}
for motifLen in $(echo ${motifLens} | sed -n 1'p' | tr ',' '\n')
do
    scriptName=${motifData}/$(basename ${motifData})_prepareSNPs_len${motifLen}.sh
    motifData_curLen=${motifData}/SNP2bed_motifLen${motifLen}
    mkdir ${motifData_curLen}
    mkdir ${motifData_curLen}/Fasta
    out=${motifData_curLen}/$(basename ${SNPQTL})_motifLen${motifLen}
    mLe=${motifLen}
    echo "source ${source_place}" > ${scriptName}
    echo "zcat ${SNPQTL} | awk '{snpend="'$1'"+1}{print \"chr\""'$2'"\"\t\""'$1'"\"\t\"snpend}' | grep -v \"chrchr\" | sort | uniq > ${out}.snpsite" >> ${scriptName}
    past_output=${out}.snpsite
    for filter in $(echo ${filters} | sed -n 1'p' | tr ',' '\n')
    do 
	echo ${filter}
	cur_output=${past_output}.filter.$(echo $(basename ${filter}) | sed 's/.bed.gz//g')
	echo "zcat ${filter} | bedtools intersect -wa -a ${past_output} -b stdin > ${cur_output}" >> ${scriptName}
	echo "rm ${past_output}" >> ${scriptName}
	past_output=${cur_output}
    done
    echo "cat ${cur_output} | awk '{mL="'$2'"-${mLe}+1}{pL="'$2'"+${mLe}}{print "'$1'"\"\t\"mL\"\t\"pL}' | grep -v \"chrc\" | sort | uniq | gzip > ${cur_output}.bed.gz" >> ${scriptName}
    echo "rm ${cur_output}" >> ${scriptName}
    chmod 711 ${scriptName}
    initialID=`qsub -l hostname=scg3-0-* -o ${scriptName}.o -e ${scriptName}.e ${scriptName} | cut -d ' ' -f 3`

    #Get fasta sequences from the individuals
    while read line
    do
	
	items=($(echo ${line} | tr ' ' '\n'))
	person=${items[0]}
	cur_script=${scriptName}_${person}_${motifLen}_script.sh
	echo "source ${source_place}" > ${cur_script}
	maternal=${items[1]}
	paternal=${items[2]}
	person_fasta=${motifData_curLen}/Fasta/${person}
	basename_regions=$(basename ${cur_output}.bed.gz) 
	echo "Processing ${person} with maternal genome ${maternal} and paternal genome ${paternal}"
	echo "cp ${cur_output}.bed.gz ${person_fasta}REGIONS.gz" >> ${cur_script}
	echo "zcat ${person_fasta}REGIONS.gz | bedtools getfasta -fi ${maternal} -bed stdin -fo ${person_fasta}.maternalSeqsIN.${basename_regions}.fa" >> ${cur_script}
	echo "zcat ${person_fasta}REGIONS.gz | bedtools getfasta -fi ${paternal} -bed stdin -fo ${person_fasta}.paternalSeqsIN.${basename_regions}.fa" >> ${cur_script}
	echo "rm ${person_fasta}REGIONS.gz" >> ${cur_script}
	echo "cat ${cur_script} >> ${cur_script}.o" >> ${cur_script}
	echo "rm ${cur_script} ${cur_script}.e ${cur_script}.o" >> ${cur_script}
	chmod 711 ${cur_script}
	qsub -hold_jid ${initialID} -l hostname=scg3-0-* -wd /srv/gsfs0/scratch/oursu/ -l h_vmem=10G -o ${cur_script}.o -e ${cur_script}.e ${cur_script}
    done < ${metadata}
done

