#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# shell script for filtlong filtering of minION reads against illumina reads
#######################################################################################
#######################################################################################
# https://github.com/rrwick/Filtlong

## installation
# ssh grace
# hpc # git is automatically available
# interactive -R "rusage[mem=4000]" -M 4000 # 4000=4GB RAM, which is the max for an hpc interactive job
# cd ~/scripts/
# git clone https://github.com/rrwick/Filtlong.git
# cd Filtlong
# module load gcc
# make -j
# PATH=$PATH:/gpfs/home/b042/scripts/Filtlong/bin
# filtlong -h

# interactive usage, but crashes a lot with only 4GB of RAM
interactive -x -R "rusage[mem=4000]" -M 4000 # 4000=4GB RAM, which is the max for an hpc interactive job
module load gcc
PATH=$PATH:/gpfs/home/b042/scripts/Filtlong/bin
PATH=$PATH:/gpfs/home/b042/scripts/seqtk
PATH=$PATH:/gpfs/home/b042/scripts/parallel-20170722/bin/
WORKDIR=~/pollen_minION/
ILMNDIR=plantskims/  # read datasets, original, not denoised, not adapter- or quality-trimmed,
NPROC=$(grep -c processor /proc/cpuinfo); echo $NPROC
filtlong -h
parallel -h
cd ${WORKDIR}
ls ${ILMNDIR}*.fq.gz
ls ${WORKDIR}minION/barcode*.fastq
MINIONREADS=10 # minION reads from a particular barcode (pollen soup)
ILMNREF=B4  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 --mean_q_weight 10 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)

zgrep -c "^" output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz


# minION 12 and Illumina B4 (343M):  5924/4 = 1481 reads
     # minION barcode 12 should have a few good matches to B4 Plantago lanceolata
# minION 12 and Illumina C12 (233M):  5924/4 = 1481 reads
     # minION barcode 12 should have a few good matches to C12 Mimulus guttatus
# minION 12 and Illumina B2:  /4 =  reads
# minION 12 and Illumina E3 (147M):  /4 =  reads
# minION 12 and Illumina B5:  0/4 =  0 reads


# minION barcode 12 should have a few good matches to B2 Ranunculus acris
# minION barcode 12 should have 0 good matches to C4 Ranunculus repens
# minION barcode 12 should have the most good matches to E3 Epilobium hirsutum
# minION barcode 12 should have 0 good matches to B5 Papaver somniferum (but is a congener of A5 so some false positives are expected)
# minION barcode 12 should have many good matches to A5 Papaver rhoeas
# minION barcode 12 should have 0 good matches to A8 Tripleurospermum maritimum

# real data from soups
#    01	02	03	04	05	06	07	08	09	10	11	12
# A2	0	0	0	0	0	0	0	0	0	1	1	0	Reseda luteola
# A3	0	0	0	0	0	0	0	0	0	1	1	0	Trifolium repens
# A5	0	0	1	1	1	0	0	0	0	1	90	2	Papaver rhoeas
# A5	1	0	0	0	0	0	0	1	9	0	0	0	Papaver rhoeas
# A6	0	0	0	0	1	0	0	0	0	1	1	0	Echium vulgare
# B2	0	0	0	0	0	0	0	0	0	1	1	1	Ranunculus acris
# B4	0	0	0	0	0	0	1	0	0	1	1	2	Plantago lanceolata
# C2	0	0	0	0	0	0	0	0	0	1	1	0	Lotus corniculatus
# C4	0	0	0	0	0	0	1	0	0	1	1	0	Ranunculus repens
# C12	0	0	0	1	0	1	0	0	0	1	1	0	Mimulus guttatus
# E3	0	1	1	0	0	1	0	9	1	90	1	5	Epilobium hirsutum
# C6	0	0	0	0	0	0	0	0	0	1	1	0	Anagallis arvensis

# loop version
########## BSUB START ###################################################
#!/bin/sh
#BSUB -q mellanox-ib  # short-eth has 24 hrs, mellanox-ib 7 days
#BSUB -n 28  # short-eth has 24 cores, add #BSUB -x # if i want an exclusive
#BSUB -J filt  # 10 chars max
#BSUB -B
#BSUB -N
#BSUB -oo filt.out
#BSUB -eo filt.err
#BSUB -R "rusage[mem=40000]" # 40000 = 40 GB RAM.
#BSUB -M 40000

. /etc/profile
module load gcc
PATH=$PATH:/gpfs/home/b042/scripts/Filtlong/bin
WORKDIR=~/pollen_minION/
ILMNDIR=plantskims/  # read datasets, original, not denoised, not adapter- or quality-trimmed,
MINIONDIR=minION/

cd ${WORKDIR}
declare -a illuminarray=(A2 A3 A5 A8 B2 B4 B5 C2 C4 C12 E3 C6)
declare -a minionarray=(10 12)

for i in "${illuminarray[@]}"
do
     for m in "${minionarray[@]}"
     do
          echo "Illumina reference $i"
          echo "minION query $j"
          filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${i}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${i}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 --mean_q_weight 10 ${WORKDIR}${MINIONDIR}barcode${m}_all_pass.fastq | gzip > output_ILMNREF_${i}_MINIONREADS_${m}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)
     done
done
########## BSUB END ###################################################

# non-loop version
########## BSUB START ###################################################
#!/bin/sh
#BSUB -q mellanox-ib  # short-eth has 24 hrs, mellanox-ib 7 days
#BSUB -n 28  # short-eth has 24 cores, add #BSUB -x # if i want an exclusive
#BSUB -J filt  # 10 chars max
#BSUB -B
#BSUB -N
#BSUB -oo filt.out
#BSUB -eo filt.err
#BSUB -R "rusage[mem=40000]" # 40000 = 40 GB RAM.
#BSUB -M 40000

. /etc/profile
module load gcc
PATH=$PATH:/gpfs/home/b042/scripts/Filtlong/bin
WORKDIR=~/pollen_minION/
ILMNDIR=plantskims/  # read datasets, original, not denoised, not adapter- or quality-trimmed,

cd ${WORKDIR}

MINIONREADS=12 # minION reads from a particular barcode (pollen soup)
ILMNREF=A2  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)


ILMNREF=A3  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)


ILMNREF=A5  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)

ILMNREF=A6  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)

ILMNREF=A8  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)

ILMNREF=B2  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)


ILMNREF=B4  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)


ILMNREF=B5  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)


ILMNREF=C2  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)


ILMNREF=C4  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)


ILMNREF=C12  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)


ILMNREF=E3  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)


ILMNREF=A5  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)

ILMNREF=C6  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_${ILMNREF}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode${MINIONREADS}_all_pass.fastq | gzip > output_ILMNREF_${ILMNREF}_MINIONREADS_${MINIONREADS}.fastq.gz  # the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)

########## BSUB END ###################################################



# parallel version # CRASHES after the first one or two Illumina libraries.  dont' know why
########## BSUB START ###################################################
#!/bin/sh
#BSUB -q mellanox-ib  # short-eth has 24 hrs, mellanox-ib 7 days
#BSUB -n 28  # short-eth has 24 cores, add #BSUB -x # if i want an exclusive
#BSUB -J filtpll  # 10 chars max
#BSUB -B
#BSUB -N
#BSUB -oo filtpll.out
#BSUB -eo filtpll.err
#BSUB -R "rusage[mem=40000]" # 40000 = 40 GB RAM.
#BSUB -M 40000

. /etc/profile
module load gcc
PATH=$PATH:/gpfs/home/b042/scripts/Filtlong/bin
PATH=$PATH:/gpfs/home/b042/scripts/parallel-20170722/bin/
WORKDIR=~/pollen_minION/
ILMNDIR=plantskims/  # read datasets, original, not denoised, not adapter- or quality-trimmed,

cd ${WORKDIR}
parallel filtlong -1 ${WORKDIR}${ILMNDIR}Plate1_{1}_*_L002_R1_val_1.fq.gz -2 ${WORKDIR}${ILMNDIR}Plate1_{1}_*_L002_R2_val_2.fq.gz --min_length 1000 --keep_percent 90 ${WORKDIR}minION/barcode{2}_all_pass.fastq ">" output_parallel_ILMNREF_{1}_MINIONREADS_{2}.fastq ::: A2 A3 A5 A6 A8 B2 B4 B5 C2 C4 C12 E3 C6 ::: 10 12  # quote the redirect because it is processed outside of parallel

parallel gzip output_parallel_ILMNREF_{1}_MINIONREADS_{2}.fastq ::: A2 A3 A5 A6 A8 B2 B4 B5 C2 C4 C12 E3 C6 ::: 10 12

# the wildcard in the illumina reference is because the index sequence (e.g. AGTTCAATC) varies, but there is only one file for each plate well reference (e.g. A1)

########## BSUB END ###################################################

# minION barcode 12 should have a few good matches to B2 Ranunculus acris
# minION barcode 12 should have 0 good matches to C4 Ranunculus repens
# minION barcode 12 should have the most good matches to E3 Epilobium hirsutum
# minION barcode 12 should have 0 good matches to B5 Papaver somniferum (but is a congener of A5 so some false positives are expected)
# minION barcode 12 should have many good matches to A5 Papaver rhoeas
# minION barcode 12 should have 0 good matches to A8 Tripleurospermum maritimum

# real data from soups
#    01	02	03	04	05	06	07	08	09	10	11	12
# A2	0	0	0	0	0	0	0	0	0	1	1	0	Reseda luteola
# A3	0	0	0	0	0	0	0	0	0	1	1	0	Trifolium repens
# A5	0	0	1	1	1	0	0	0	0	1	90	2	Papaver rhoeas
# A5	1	0	0	0	0	0	0	1	9	0	0	0	Papaver rhoeas
# A6	0	0	0	0	1	0	0	0	0	1	1	0	Echium vulgare
# B2	0	0	0	0	0	0	0	0	0	1	1	1	Ranunculus acris
# B4	0	0	0	0	0	0	1	0	0	1	1	2	Plantago lanceolata
# C2	0	0	0	0	0	0	0	0	0	1	1	0	Lotus corniculatus
# C4	0	0	0	0	0	0	1	0	0	1	1	0	Ranunculus repens
# C12	0	0	0	1	0	1	0	0	0	1	1	0	Mimulus guttatus
# E3	0	1	1	0	0	1	0	9	1	90	1	5	Epilobium hirsutum
# C6	0	0	0	0	0	0	0	0	0	1	1	0	Anagallis arvensis
