

samples=(`zcat /srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/results/Pickrell/GSE19480_series_matrix.txt.gz | sed -n '30p' | sed 's/"//g'`)
files=(`zcat /srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/results/Pickrell/GSE19480_series_matrix.txt.gz | sed -n '63p' | sed 's/"//g'`)
fastqdumpcode=/srv/gs1/software/sratoolkit/sratoolkit.2.3.4-2-ubuntu64/bin/fastq-dump
pickrelldir=/srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/results/Pickrell/
pickrelldir=/srv/gsfs0/projects/kundaje/users/oursu/histoneQTLproject/RNAseq/results/Pickrell/
for ((i=1; i<${#samples[*]}; i++)); 
do
    echo "here"
    ind=$(echo ${samples[$i]} | sed 's/_2//g' | sed 's/_yale//g' | sed 's/_argonne//g')
    outlab=$(echo ${samples[$i]} | sed 's/_yale/.yale/g' | sed 's/_argonne/.argonne/g')
    outname=${pickrelldir}/$ind/${outlab}.fastq.gz
    if [[ ${outname} == *argonne* ]]
	then
	echo "Name Found------------------------------"
	outscript=${outname}_script2.sh
	mkdir ${pickrelldir}/$ind
	echo $ind
	echo ${samples[$i]}
	echo ${files[$i]}
	echo "cd ${pickrelldir}/${ind}/" > ${outscript}
	#echo "wget -P ${samples[$i]} -r -nd ${files[$i]}" >> ${outscript}
        echo "${fastqdumpcode} --gzip -Z -F ${samples[$i]}/*.sra > ${samples[$i]}.fastq.gz" >> ${outscript}
	echo ${outscript}
	chmod 711 ${outscript}
	qsub -l h_vmem=3G -o ${outscript}.o -e ${outscript}.e ${outscript}
    fi
done
