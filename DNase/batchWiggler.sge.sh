#!/bin/bash
# ###########################################
# AUXILIARY FUNCTIONS
# ###########################################

# Checks if a file or directory exists and exits if it fails
FUNCdoesFileExist()
{
    local VarName=$1
    local fileName=$2    
    if [[ ! -e "${fileName}" ]]
    then
        echo "ERROR in variable ${VarName}:${fileName} - File/Directory does not exist"
        exit 1
    fi
}

# ###########################################
# MAIN
# ###########################################

if [[ "$#" -lt 13 ]]
then
	echo "$(basename $0): Will create signal files from a set of tagAlign/BAM files and corresponding extension lengths" >&2
	echo "USAGE: $(basename $0) <PairFile> <FraglenFile> <inputDir> <outputDir> <replaceFlag> <memory> <oformat> <kernelName> <normFlag> <mapFilter> <seqDir> <umapDir>" >&2
	echo ' <PairFile>: Col1: replicate alignment file names separated by ; Col2: OutputFilePrefix' >&2
	echo ' <FraglenFile>: tab delimited file where Col1: FileName and Col2: fragment length (2*tag-shift) and Col3: smoothing window size' >&2
	echo ' <inputDir>: Root directory under which all the align files are located' >&2
	echo ' <outputDir>: directory where you want the output files to be' >&2
	echo ' <replaceFlag>: TRUE/FALSE If set to TRUE then files will be overwritten else skipped' >&2
	echo ' <memory>: in GB (typically 4 or 6 works fine)' >&2
	echo ' <oformat>: mat,bg,wig (Preferably use bg. It is the most efficient to convert to bigwig)' >&2
	echo ' <kernelName>: rectangular,triangular,epanechnikov,biweight,triweight,cosine,gaussian,tukey (For most data type use epanechnikov. For MNase-seq use triweight. If you just want tag extension use rectangular)' >&2
	echo ' <normFlag>: 0: counts, 5: foldchange' >&2
	echo ' <mapFilter>: 0 allows all locations' >&2
	echo ' <seqDir>: directory containing chromosome sequences e.g. /srv/gsfs0/projects/kundaje/users/akundaje/projects/encode/data/byDataType/sequence/encodeHg19Male' >&2
	echo ' <umapDir>: directory containing mappability files e.g. /srv/gsfs0/projects/kundaje/users/akundaje/projects/encode/data/byDataType/umap/hg19_allmappable/globalmap_k1tok1000' >&2
	echo ' <chrSizes>' >&2
	echo ' Dependencies: Make sure the align2rawsignal binary is in your $PATH' >&2
	echo ' Dependencies: $TMP should point to a temporary directory' >&2
	echo ' Dependencies: The code current points to Anshuls MCR installation on the SCG3 cluster. If you are running this on a different machine have the MCR installed on that machine and change $MCRROOT in the code below' >&2
	exit 1
fi

PAIRFILE=$1 # File containing ChIP file names and output file prefix
FUNCdoesFileExist 'PAIRFILE' "${PAIRFILE}"

FRAGLENFILE=$2 # File containing fragment lengths and smoothing sizes
FUNCdoesFileExist 'FRAGLENFILE' "${FRAGLENFILE}"

IDIR=$3 # Input root directory containing tagAlign/BAM files
FUNCdoesFileExist 'IDIR' "${IDIR}"
printf "ALIGNDIR: %s\n" ${IDIR}

ODIR=$4 # Output directory
FUNCdoesFileExist 'ODIR' "${ODIR}"
printf "OUTPUTDIR: %s\n" ${ODIR}

REPLACEFLAG=$5 # If set to TRUE then if output file exists, it is replaced. Else it is skipped
printf "REPLACEFLAG: %s\n\n" ${REPLACEFLAG}

MEM=$6

OFORMAT=$7

MAXJOBS=30

KERNELNAME=$8

JOBGROUP="/wiggler${RANDOM}"
MAXJOBS=30
#bgadd -L ${MAXJOBS} $JOBGROUP

NORMFLAG=$9 # Normalization flag
printf "NORMFLAG: %d\n" ${NORMFLAG}

MAPFILTER=${10}
printf "MAPFILTER: %s\n" ${MAPFILTER}

printf "OUTPUTFORMAT: %s\n" ${OFORMAT}
if [[ "${OFORMAT}" == 'bg' ]]
then
	OFEXT=".norm${NORMFLAG}.rawsignal.bedgraph" # Output file extension
else
	OFEXT=".norm${NORMFLAG}.rawsignal.${OFORMAT}" # Output file extension
fi
LOGEXT=".norm${NORMFLAG}.rawsignal.log" # log file extension

CHR_DIR=${11}
FUNCdoesFileExist 'CHR_DIR' "${CHR_DIR}"

UNIQ_DIR=${12}
FUNCdoesFileExist 'UNIQ_DIR' "${UNIQ_DIR}"

CHRSIZE=${13}

BASECMD="align2rawsignal -mm=${MEM} -of=${OFORMAT} -n=${NORMFLAG} -f=${MAPFILTER} -k=${KERNELNAME} -s=${CHR_DIR} -u=${UNIQ_DIR}" # Base ta2rs command                                                                                                                        

count=0                                                                                             
while read line 
do
    chipfiles=$(echo ${line} | awk '{print $1}')
    outputPrefix=$(echo ${line} | awk '{print $2}')

    # ====================================	    
    # Generate logfile name
    # ====================================
    LOGFILE="${ODIR}/${outputPrefix}${LOGEXT}"    
    OFNAME="${ODIR}/${outputPrefix}${OFEXT}"
    CMD="$BASECMD -v=${LOGFILE}" # update command with --v
    
    # ====================================   
    # Get files and fragment lengths corresponding to the root
    # ====================================
    nf=0
    for i in $(echo ${chipfiles} | sed -r 's/;/\n/g')
    do
	c=$((c+1))
	filepath=$(find ${IDIR} -name $i | head -1)
	if [[ ! -e ${filepath} ]]
	then
	    echo "WARNING: File $i not found" >&2
	    nf=1
	    continue
	else
	    fragline=$(grep -F $i ${FRAGLENFILE} | head -1)
	    if [[ -z ${fragline} ]]
	    then
		echo "WARNING: No fragment length information for file $i" >&2
		nf=1
	    else
		fraglen=$(echo ${fragline} | awk '{print $2}')
		smoothlen=$(echo ${fragline} | awk '{print $3}')
		CMD="$CMD -i=${filepath} -l=${fraglen}"
	    fi
	fi
    done
    CMD="$CMD -w=${smoothlen}"
    if [[ ${nf} -eq 1 ]]
    then
        echo "Error for ${line}" >&2
	continue
    fi
    
    # ====================================
    # Generate outputfile name
    # ====================================
    # Update command with --ofile and optionally --output-max-tags
    if [[ "${OFORMAT}" == 'mat' ]]
    then
        CMD="$CMD -o=${OFNAME}" # Update command with --ofile
	#MAXTAGSFNAME=$(echo ${OFNAME} | sed -r 's/\.rawsignal\./.rawsignal.maxtags./g')
    	#CMD="$CMD --output-max-tags=${MAXTAGSFNAME}"
    else
    	OFNAME="${OFNAME}.gz"
    	CMD="$CMD | grep -v Warning | gzip -c 1> ${OFNAME}"
    fi

    # ====================================
    # Check if output file exists and if REPLACEFLAG is not set then skip
    # ====================================
    if [[ -e "${OFNAME}" && "${REPLACEFLAG}" != 'TRUE' ]]
    then
        echo "Output File Exists: Skipping ${OFNAME}"
        printf "%s\t%s\tSkipped\n" $(date +%D) ${OFNAME} >> skippedDatasets.log
        continue
    fi

    export PATH=/srv/gs1/software/samtools/samtools-0.1.19/bin/:$PATH
    # Create temp submit/run script
    TMP=${ODIR}/${outputPrefix}_TMP2
    mkdir ${TMP}
    SUBMITFILE="${ODIR}/${outputPrefix}.sh"
    echo '#!/bin/bash' > ${SUBMITFILE}
    #echo "export PATH=/srv/gs1/software/samtools/samtools-0.1.19/bin/:$PATH" >> ${SUBMITFILE}
    #echo 'export MCRROOT=/srv/gsfs0/projects/kundaje/users/akundaje/local/lib/mcr/2010b/v714' >> ${SUBMITFILE}
    #echo 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/runtime/glnxa64' >> ${SUBMITFILE}
    #echo 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64' >> ${SUBMITFILE}
    #echo 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64' >> ${SUBMITFILE}
    #echo 'MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64' >> ${SUBMITFILE}
    #echo 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads' >> ${SUBMITFILE}
    #echo 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server' >> ${SUBMITFILE}
    #echo 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}' >> ${SUBMITFILE}
    #echo 'XAPPLRESDIR=${MCRROOT}/X11/app-defaults' >> ${SUBMITFILE}
    #echo 'export LD_LIBRARY_PATH' >> ${SUBMITFILE}
    #echo 'export XAPPLRESDIR' >> ${SUBMITFILE}
    #echo "source /srv/gsfs0/projects/kundaje/commonRepository/src/lab_bashrc" >> ${SUBMITFILE}


    TMP=/srv/gsfs0/scratch/oursu
    mkdir ${TMP}
    echo "#source  /srv/gsfs0/projects/kundaje/commonRepository/src/lab_bashrc" >> ${SUBMITFILE}
    echo "#activateMCR" >> ${SUBMITFILE}
    echo "#TEMPDIR=${TMP}/"'temp_${RANDOM}${RANDOM}' >> ${SUBMITFILE} # Create temp directory for source files
    echo '#mkdir $TEMPDIR' >> ${SUBMITFILE}
    echo '#export MCR_CACHE_ROOT=${TEMPDIR}' >> ${SUBMITFILE}
    echo "#cp $(which align2rawsignal) "'${TEMPDIR}/' >> ${SUBMITFILE}
    echo '#cd ${TEMPDIR}' >> ${SUBMITFILE}
    echo "#./${CMD}" >> ${SUBMITFILE}
    echo "module load bedtools/2.17.0" >> ${SUBMITFILE} 
    echo "zcat ${ODIR}/${outputPrefix}.norm5.rawsignal.bedgraph.gz | slopBed -i stdin -g ${CHRSIZE} -b 0 | /home/oursu/devtools/bedClip stdin ${CHRSIZE} ${TMP}/${outputPrefix}.norm5.rawsignal.slopclip.bedgraph" >> ${SUBMITFILE} 
    echo "/home/oursu/devtools/toolsForBedGraphics/bedGraphToBigWig ${TMP}/${outputPrefix}.norm5.rawsignal.slopclip.bedgraph ${CHRSIZE} ${ODIR}/${outputPrefix}.norm5.rawsignal.bigwig" >> ${SUBMITFILE}
    #echo "rm ${ODIR}/${outputPrefix}.norm5.rawsignal.slopclip.bedgraph ${ODIR}/${outputPrefix}.norm5.rawsignal.bedgraph.gz" >> ${SUBMITFILE}
    echo 'cd ..' >> ${SUBMITFILE}
    echo '#rm -rf ${TEMPDIR}' >> ${SUBMITFILE}
    chmod 755 ${SUBMITFILE}
    echo ${SUBMITFILE}
    # Submit/Run the script
    #qsub -cwd -q standard -w e -V -j y -o ${LOGFILE} -e "${LOGFILE}.err" -l h_vmem="${MEM}G" -l h_rt=6:00:00 -N ${outputPrefix} "${SUBMITFILE}"
    qsub -cwd -q standard -w e -V -j y -l h_vmem=5G -o ${LOGFILE} -e "${LOGFILE}.err" -N ${outputPrefix} "${SUBMITFILE}"
    #rm ${SUBMITFILE}

done < ${PAIRFILE}
