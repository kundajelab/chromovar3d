cmd='qsub -l h_vmem=10G -l h_rt=20:00:00 -o /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/src/download_mapped_reads.sh.o -e /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/src/download_mapped_reads.sh.e /srv/gs1/projects/snyder/jzaugg/histoneQTL/DNase/src/download_mapped_reads.sh'
echo $cmd
$cmd