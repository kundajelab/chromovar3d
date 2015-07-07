#!/bin/bash
# ========================================
# Read in arguments and check for errors
# ========================================
if [[ "$#" -lt 3 ]]
    then
    echo "Usage: $(basename $0) <iDir> <pairFile> <oDir> <OPT:genome> <OPT:chrSizes> <OPT:memory> <OPT:preComputedFragLenFile>" >&2
    echo '<iDir>: input directory whose subdirectories contain mapped ChIP and control data' >&2
    echo '<pairFile>: input file containing pairs of ChIP (column 1) and control dataset names (column 2).' >&2 
    echo '            Multiple ChIP and control files can be present seperated by ;. They will be merged before passing them to MACS' >&2
    echo '<oDir>: output directory' >&2
    echo '<genome>: (OPTIONAL) genome size. Default: hs' >&2
    echo '<chrSizes>: (OPTIONAL) tag delimited chrName\tsize: Default : $GENOMESIZEDIR/hg19.genome' >&2
    echo '<memory>: (OPTIONAL) memory limit 0 (means default), if non-zero then it represents X GB' >&2
    echo '<preComputedFragLenFile>: (OPTIONAL) File containing estimated fragment length for each dataset. Col1: dataset name [tab] Col2:frag length' >&2		
    exit 1
fi

IDIR=$1
if [[ ! -d "${IDIR}" ]]; then echo "ERROR: Argument <iDir> ${IDIR} does not exist" >&2 ; exit 1; fi

IFILE=$2
if [[ ! -f "${IFILE}" ]]; then echo "ERROR: Argument <pairFile> ${IFILE} does not exist" >&2 ; exit 1; fi

ODIR=$3
if [[ ! -d "${ODIR}" ]]; then echo "ERROR: Argument <oDir> ${ODIR} does not exist" >&2 ; exit 1; fi

GENOMESIZE='hs'
if [[ "$#" -ge 4 ]]; then GENOMESIZE=$4 ; fi

CHRSIZE="${GENOMESIZEDIR}/hg19.genome"
if [[ "$#" -ge 5 ]]; then CHRSIZE=$5 ; fi
if [[ ! -e "${CHRSIZE}" ]]; then echo "ERROR: Argument <chrSizes> ${CHRSIZE} does not exist" >&2 ; exit 1; fi

MEM=0
if [[ "$#" -ge 6 ]]; then MEM=$6 ; fi
#memLimit=$(( MEM * 1024 ))
memLimit=${MEM}

FRAGLENFILE='NULL'
if [[ "$#" -ge 7 ]]; then FRAGLENFILE=$7 ; fi
if [[ ${FRAGLENFILE} != 'NULL' && ! -f "${FRAGLENFILE}" ]]
then
    echo "ERROR: Argument <iSizeFile> ${FRAGLENFILE} does not exist" >&2
    exit 1
fi

# Check if MACS is in the path
#if [[ -z $(which macs2) ]]; then echo 'ERROR: MACS executable not in $PATH' >&2; exit 1; fi

# Maximum number of jobs to run at a time
#JOBGRPID="/compbio/long/macs2signal${RANDOM}"

# ========================================
# Read pairFile line by line
# Create shell script
# Submit to cluster
# ========================================


while read inputline
do

  # -------------------------
  # Operate on ChIP files
  # -------------------------
  chip=$(echo "$inputline" | awk 'BEGIN{FS="\t| "}{print $1}') # extract chip file
  chipstub=$(echo ${chip} | sed -r -e 's/([^;]*);.*$/\1/g' -e 's:^.*/::g' -e 's/\.gz$//g') # Get first replicate name, remove path information and remove .gz extension
  [[ -n $(echo ${chip} | grep ';') ]] && chipstub=$(echo ${chipstub} | sed -r 's/_[0-9]_+/_0_/g') # If multiple file names then replace the Rep[\d]+ in chipstub with Rep0

  # --------------------------
  # Operate on controls
  # --------------------------
  noControl=0
  control=$(echo "$inputline" | awk 'BEGIN{FS="\t| "}{print $2}') # extract control file
  [[ -z ${control} ]] && noControl=1
  if [[ ${noControl} == '0' ]]
      then
      controlstub=$(echo ${control} | sed -r -e 's/([^;]*);.*$/\1/g' -e 's:^.*/::g' -e 's/\.gz$//g') # Get first replicate name, remove path information and remove .gz extension
      [[ -n $(echo ${control} | grep ';') ]] && controlstub=$(echo "${controlstub}" | sed -r 's/_[0-9]_+/_0_/g') # If multiple file names then replace the Rep[\d]+ in chipstub with Rep0
  else
      controlstub='NoControl'
      echo ${controlstub}
  fi
  echo 'then'
  echo ${noControl}
  echo 'control'
  echo ${controlstub}

  # ----------------------
  # Set Output file names
  # ----------------------
  outFile="${ODIR}/${chipstub}_VS_${controlstub}.out"
  errFile="${ODIR}/${chipstub}_VS_${controlstub}.err"
  peakFile="${ODIR}/${chipstub}_VS_${controlstub}"

  #if [[ -e "${outFile}" || -e "${errFile}" || -e "${peakFile}_peaks.narrowPeak" ]]
  #    then
  #    echo "Skipping ${peakFile}" >&2
  #    continue
  #fi

  chipstub=$(echo ${chipstub} | sed -r 's/\.bam$/.tagAlign/g')  # Change ChIPstub and controlstub to tagAlign extensions
  chip=$(echo "${chip}" | sed -r 's/\.bam(\.gz)?(;|$)/.bam\1.unique.tagAlign.gz\2/g') # Change file extensions from bam|bam.gz to bam.unique.tagAlign.gz|bam.gz.unique.tagAlign.gz

  if [[ ${noControl} == '0' ]]
      then
      controlstub=$(echo ${controlstub} | sed -r 's/\.bam$/.tagAlign/g') # Change controlstub to tagAlign extensions
      control=$(echo "${control}" | sed -r 's/\.bam(\.gz)?(;|$)/.bam\1.unique.tagAlign.gz\2/g')  # Change file extensions from bam|bam.gz to bam.unique.tagAlign.gz|bam.gz.unique.tagAlign.gz
  fi

  # -------------------------
  # Search for full file path in IDIR
  # ; are replaced by <space>
  # -------------------------
  nchip=$(echo ${chip} | sed -r 's/;/\n/g' | wc -l) # number of chip files
  chip=$(echo "${chip}" | sed -r 's/;/\n/g' | xargs -I fname find "${IDIR}" -name fname -printf "%p ")
  nchipfound=$(echo ${chip} | sed -r 's/ /\n/g' | wc -l) # number of chip files that were found in IDIR
  if [[ ${nchipfound} -ne ${nchip} ]]; then echo "ERROR: Some of the ChIP files ${chip} were not found in $IDIR" >&2 ; continue ; fi
  
  if [[ ${noControl} == '0' ]]
      then
      ncontrol=$(echo ${control} | sed -r 's/;/\n/g' | wc -l) # number of control files
      control=$(echo "${control}" | sed -r 's/;/\n/g' | xargs -I fname find "${IDIR}" -name fname -printf "%p ")
      ncontrolfound=$(echo ${control} | sed -r 's/ /\n/g' | wc -l) # number of control files that were found in IDIR
      if [[ ${ncontrolfound} -ne ${ncontrol} ]]; then echo "ERROR: Some of the control files ${control} were not found in $IDIR" >&2 ; continue ; fi
  fi

  # --------------------------
  # Get fragLens corresponding to ChIP files
  # --------------------------
  fraglen=0
  fcount=0
  if [[ ${FRAGLENFILE} != 'NULL' ]]
  then
      echo "chip"
      echo ${chip}
      for currFile in $(echo ${chip} | sed -r 's/;/\n/g')
      do
	  echo "current file"
	  echo ${currFile}
	  currBase=$(echo $(basename ${currFile} | sed -r 's/(\.bam.*$)|(\.tagAlign.*$)/\\./g'))
	  if grep -q ${currBase} ${FRAGLENFILE}
	  then
	      currFragLen=$(grep ${currBase} ${FRAGLENFILE} | awk '{print $2}')
	      fraglen=$((fraglen + currFragLen))
	      fcount=$((fcount + 1))
	  else
	      echo "WARNING: No fragment length found corresponding to file ${currFile}" >&2
	  fi
	  echo "frag len"
	  echo ${fraglen}
      done
      [[ ${fcount} -eq 0 ]] && fcount=1
      fraglen=$((fraglen / (2 * fcount)))
      if [[ ${fraglen} -eq 0 ]]
      then
	  echo "ERROR: Fragment length is 0 due to missing file names in ${FRAGLENFILE}" >&2
	  continue
      fi
  fi
  
  # -------------------------
  # Initialize submit script
  # -------------------------
  #if [[ $(bjobs -g "${JOBGRPID}" 2> /dev/null | wc -l) -gt 60 ]]
  #    then
  #    sleep 30s
  #fi

  #scriptName="temp${RANDOM}${RANDOM}.sh" # script name
  scriptName="${peakFile}_MACS2_script.sh"
  echo '#!/bin/bash' > ${scriptName}
  echo 'module load MACS2/2.0.10' >> ${scriptName}
  echo 'tmpdir="${TMP}/tmp${RANDOM}_${RANDOM}"' >> ${scriptName}
  echo 'mkdir ${tmpdir}' >> ${scriptName}
  
  # -------------------------
  # Create temp copies of ChIP and control
  # gunzip and concatenate multiple files if necessary
  # -------------------------	
  echo 'combchip="${tmpdir}/'"${chipstub}"'"' >> ${scriptName}
  #echo 'if [[ -f "${combchip}" ]]; then rm -rf "${combchip}"; fi' >> ${scriptName}
  echo "echo Combining ChIP replicates: ${chip}" >> ${scriptName}
  echo "gunzip -c ${chip} >> "'"${combchip}"' >> ${scriptName}
  echo 'chipReads=$(wc -l ${combchip} | awk '"'"'{printf "%f", $1/1000000}'"')" >> ${scriptName}
  echo 'sval=${chipReads}' >> ${scriptName}
  
  if [[ ${noControl} == '0' ]]
      then
      echo 'combcontrol="${tmpdir}/'"${controlstub}"'"' >> ${scriptName}
      #echo 'if [[ -f "${combcontrol}" ]]; then rm -rf "${combcontrol}"; fi' >> ${scriptName}
      echo "echo Combining Control replicates: ${control}" >> ${scriptName}
      echo "gunzip -c ${control} >> "'"${combcontrol}"' >> ${scriptName}
      echo 'controlReads=$(wc -l ${combcontrol} | awk '"'"'{printf "%f", $1/1000000}'"')" >> ${scriptName}
      echo 'sval=$(echo "${chipReads} ${controlReads}" | awk '"'"'$1>$2{printf "%f",$2} $1<=$2{printf "%f",$1}'"')" >> ${scriptName}
  fi

  # -------------------------
  # Complete script
  # -------------------------
  if [[ ${noControl} == '0' ]]
      then
      if [[ ${FRAGLENFILE} == 'NULL' || ${fraglen} -eq 0 ]]
      then
          echo 'macs2 callpeak -t "${combchip}" -c "${combcontrol}" -f BED'" -n ${peakFile} -g ${GENOMESIZE} -p 1e-2 --nomodel --shiftsize 73 -B --SPMR" >> ${scriptName}
      else
	  echo 'macs2 callpeak -t "${combchip}" -c "${combcontrol}" -f BED'" -n ${peakFile} -g ${GENOMESIZE} -p 1e-2 --nomodel --shiftsize ${fraglen} -B --SPMR" >> ${scriptName}
      fi
  else
      #continue
      #call without control (for DNase)
      echo 'macs2 callpeak -t "${combchip}" -f BED'" -n ${peakFile} -g ${GENOMESIZE} -p 1e-2 --nomodel --shiftsize ${fraglen} -B --SPMR" >> ${scriptName}
  fi
  
  echo 'rm -rf ${tmpdir}' >> ${scriptName}
  echo "rm -f ${peakFile}_peaks.xls ${peakFile}_peaks.bed ${peakFile}_summits.bed" >> ${scriptName}  

  # foldchange bedgraph
  if [[ ! -e "${peakFile}.fc.signal.bigwig" ]]
  then
      #if [[ ${noControl} == '0' ]]
      #then
      echo "macs2 bdgcmp -t ${peakFile}_treat_pileup.bdg -c ${peakFile}_control_lambda.bdg --o-prefix ${peakFile} -m FE" >> ${scriptName}
      echo "/srv/gs1/software/bedtools/bedtools-2.17.0/bin/slopBed -i ${peakFile}_FE.bdg -g ${CHRSIZE} -b 0 | /home/oursu/devtools/bedClip stdin ${CHRSIZE} ${peakFile}.fc.signal.bedgraph" >> ${scriptName}
      echo "rm -f ${peakFile}_FE.bdg" >> ${scriptName}
      echo "/home/oursu/devtools/toolsForBedGraphics/bedGraphToBigWig ${peakFile}.fc.signal.bedgraph ${CHRSIZE} ${peakFile}.fc.signal.bigwig" >> ${scriptName}
      echo "rm -f ${peakFile}.fc.signal.bedgraph" >> ${scriptName}
      #else
	  #make a bigwig for the dnase - TODO - make more similar to above
	  #echo "/srv/gs1/software/bedtools/bedtools-2.17.0/bin/slopBed -i ${peakFile}_treat_pileup.bdg -g ${CHRSIZE} -b 0 | /home/oursu/devtools/bedClip stdin ${CHRSIZE} ${peakFile}.bedgraph" >> ${scriptName}
	  #echo "/home/oursu/devtools/toolsForBedGraphics/bedGraphToBigWig ${peakFile}.bedgraph ${CHRSIZE} ${peakFile}.bigwig" >> ${scriptName}
	  #echo "macs2 bdgcmp -t ${peakFile}_treat_pileup.bdg -o-prefix ${peakFile} -m FE" >> ${scriptName}
          #echo "/srv/gs1/software/bedtools/bedtools-2.17.0/bin/slopBed -i ${peakFile}_FE.bdg -g ${CHRSIZE} -b 0 | bedClip stdin ${CHRSIZE} ${peakFile}.fc.signal.bedgraph" >> ${scriptName}
          #echo "rm -f ${peakFile}_FE.bdg" >> ${scriptName}
          #echo "/home/oursu/devtools/toolsForBedGraphics/bedGraphToBigWig ${peakFile}.fc.signal.bedgraph ${CHRSIZE} ${peakFile}.fc.signal.bigwig" >> ${scriptName}
          #echo "rm -f ${peakFile}.fc.signal.bedgraph" >> ${scriptName}
      #fi
  fi

  # PVAL bedgraph
  if [[ ! -e "${peakFile}.pval.signal.bigwig" ]]
  then
      #if [[ ${noControl} == '0' ]]
      #	 then
      echo "macs2 bdgcmp -t ${peakFile}_treat_pileup.bdg -c ${peakFile}_control_lambda.bdg --o-prefix ${peakFile} -m ppois -S "'${sval}' >> ${scriptName}
      echo "/srv/gs1/software/bedtools/bedtools-2.17.0/bin/slopBed -i ${peakFile}_ppois.bdg -g ${CHRSIZE} -b 0 | /home/oursu/devtools/bedClip stdin ${CHRSIZE} ${peakFile}.pval.signal.bedgraph" >> ${scriptName}
      echo "rm -rf ${peakFile}_ppois.bdg" >> ${scriptName}
      echo "/home/oursu/devtools/toolsForBedGraphics/bedGraphToBigWig ${peakFile}.pval.signal.bedgraph ${CHRSIZE} ${peakFile}.pval.signal.bigwig" >> ${scriptName}
      echo "rm -f ${peakFile}.pval.signal.bedgraph" >> ${scriptName}     
      #fi
  fi
  #echo "rm -f ${peakFile}_treat_pileup.bdg ${peakFile}_control_lambda.bdg" >> ${scriptName}
  #echo "rm -f ${peakFile}.bedgraph'    
  # -------------------------
  # Submit script
  # -------------------------
  chmod 755 ${scriptName}
  cp ${scriptName} ${outFile}
  echo '======================================================================' >> ${outFile}
  qsub -l mem_free=${MEM}G -l h_vmem=${MEM}G -l h_rt=20:00:00 -o ${outFile} -e ${errFile} ${scriptName} 

  #if [[ ${MEM} -eq 0 ]]
  #then
  #    bsub -P compbiofolk -q compbio-week -g "${JOBGRPID}" -J "${chipstub}" -o ${outFile} -e ${errFile} < ${scriptName}
  #else
  #    bsub -P compbiofolk -q compbio-week -g "${JOBGRPID}" -J "${chipstub}" -M "${memLimit}" -R "rusage[mem=${memLimit}]" -o ${outFile} -e ${errFile} < ${scriptName}
  #fi
  
  # -------------------------
  # Delete temporary script
  # -------------------------	
  #rm "${scriptName}"
  #I will keep the script.
  
done < "${IFILE}"

