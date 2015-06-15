module load bedtools/2.19.1
module load r/3.0.1
export CODEDIR=/srv/gsfs0/projects/kundaje/users/oursu/code/chromovar3d/MotifAnalysis/
export PATH=${CODEDIR}:$PATH
export TESTDATA=${CODEDIR}testFiles
export HOMERDIR=/srv/gsfs0/projects/kundaje/users/oursu/code/motifs/HOMER/bin/
export R_WITH_PACKAGES=/home/oursu/devtools/R-3.0.2/bin/Rscript
enrichRcode=/srv/gsfs0/projects/kundaje/users/oursu/code/genome_utils/Features/TFs/overlapEnrichment/computeOverlapEnrichment_significance.R