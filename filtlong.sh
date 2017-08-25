#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# shell script for filtlong filtering of minION reads against illumina reads
#######################################################################################
#######################################################################################

# installation
ssh grace
hpc # git is automatically available
interactive -R "rusage[mem=4000]" -M 4000 # 4000=4GB RAM, which is the max for an hpc interactive job
cd ~/scripts/
git clone https://github.com/rrwick/Filtlong.git
cd Filtlong
module load gcc
make -j
PATH=$PATH:/gpfs/home/b042/scripts/Filtlong/bin
filtlong -h

# interactive usage
interactive -R "rusage[mem=4000]" -M 4000 # 4000=4GB RAM, which is the max for an hpc interactive job
module load gcc
PATH=$PATH:/gpfs/home/b042/scripts/Filtlong/bin
PATH=$PATH:/gpfs/home/b042/scripts/seqtk
PATH=$PATH:/gpfs/home/b042/scripts/parallel-20170722/bin/parallel
WORKDIR=~/pollen_minION/
ILMNDIR=plantskims/prepare_forBWA/  # read datasets, bfc denoised, adapter trimmed, R1 and R2 combined,
NPROC=$(grep -c processor /proc/cpuinfo); echo $NPROC
filtlong -h
cd ${WORKDIR}
ls ${ILMNDIR}*.fq.gz
MINIONREADS=12 # minION reads from a particular barcode (pollen soup)
ILMNREF=A5  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

# minION barcode 12 should have many good matches to A5 Papaver rhoeas
# minION barcode 12 should have 0 good matches to B5 Papaver somniferum (but is a congener of A5 so some false positives are expected)
# minION barcode 12 should have 0 good matches to A8 Tripleurospermum maritimum
# minION barcode 12 should have the most good matches to E3 Epilobium hirsutum

# example code, can use just one illumina file
# filtlong -1 illumina_1.fastq.gz -2 illumina_2.fastq.gz --min_length 1000 --keep_percent 90 --target_bases 500000000 --trim --split 250 input.fastq.gz | gzip > output.fastq.gz

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_CCTCAGAGA_bfc_trimed.fq.gz --min_length 1000 --keep_percent 90 --target_bases 500000000 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz

zgrep -c "^" output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz
