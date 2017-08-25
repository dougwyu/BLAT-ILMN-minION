#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# shell script for BLAT mapping illumina reads to minION reads
#######################################################################################
#######################################################################################

# To do:
# 1. parallel:  run all 50 skims against all 12 barcodes
# 2. pipe filter psl file to accept only those with > 220 matches
# 3. try filtlong as another way to filter minION reads by matches to illumina files: https://github.com/rrwick/Filtlong

# idea taken from NaS pipeline from Madoui et al. (2015, BMC Genomics)
# https://github.com/institut-de-genomique/NaS
# https://doi.org/10.1186/s12864-015-1519-z
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1519-z

# example commands taken from the NaS pipeline ()
              # blat -tileSize=$TILE -stepSize=$STEP -noHead barcode12 Plate1_A5_CCTCAGAGA_L002 output.psl
              # cat $OUTPUT_DIR/tmp/ILMN_reads.fa | parallel -j $NB_PROC --cat --pipe --block 10M --recstart ">" "$BLAT -tileSize=$TILE -stepSize=$STEP -noHead $NANO_READS {} $OUTPUT_DIR/tmp/psl/blat-alignment.job{#}.tile$TILE.step$STEP.psl" 2>$OUTPUT_DIR/tmp/blat-alignment.stderr
                            # This will BLAT map reads.  $TILE=10 and $STEP=5
              # mkdir output_dir
              # cat your_psl_file.psl | awk -v PFX=output_dir '{ file=PFX"/"$14".psl"; print $0>file; }'
                            # This will create one file per MinION read with the name of Illumina reads that mapped to this particular read.

# PSL file output
# http://www.ensembl.org/info/website/upload/psl.html
# Fields are space-separated, and all 21 are required.
#
              # matches - Number of matching bases that aren't repeats.
              # misMatches - Number of bases that don't match.
              # repMatches - Number of matching bases that are part of repeats.
              # nCount - Number of 'N' bases.
              # qNumInsert - Number of inserts in query.
              # qBaseInsert - Number of bases inserted into query.
              # tNumInsert - Number of inserts in target.
              # tBaseInsert - Number of bases inserted into target.
              # strand - defined as + (forward) or - (reverse) for query strand. In mouse, a second '+' or '-' indecates genomic strand.
              # qName - Query sequence name.
              # qSize - Query sequence size.
              # qStart - Alignment start position in query.
              # qEnd - Alignment end position in query.
              # tName - Target sequence name.
              # tSize - Target sequence size.
              # tStart - Alignment start position in target.
              # tEnd - Alignment end position in target.
              # blockCount - Number of blocks in the alignment.
              # blockSizes - Comma-separated list of sizes of each block.
              # qStarts - Comma-separated list of start position of each block in query.
              # tStarts - Comma-separated list of start position of each block in target.


#### Interactive version
ssh grace
ssh hpc.uea.ac.uk
# interactive -x # exclusive job, which is alternative to interactive -R
interactive -R "rusage[mem=4000]" -M 4000 # 4000=4GB RAM, which is the max for an hpc interactive job
# module load blat # v. 36
PATH=$PATH:/gpfs/home/b042/scripts/pblat/ # pblat v.35
WORKDIR=~/pollen_minION/
MINIONDIR=minION/
ILMNFAS=plantskims/prepare_forBWA/prepare_forblast/
BLATDIR=blatout/
TILE=10
STEP=5
MATCH=213 # 85% of 251 is 213.  in other words, I am asking for the illumina read to have at least 85% similarity to the minION read
NPROC=$(grep -c processor /proc/cpuinfo)

cd ${WORKDIR}
declare -a illuminarray=(A2 A3 A5 B2 B4 C2 C4 C6 C12 E3) # # Illumina reads = plant skims that have been adapter-trimmed, bfc-denoised, merged, dereplicated, and converted to fasta
declare -a minionarray=(11 12) # minION reads, all barcodes that successfully sequenced

# e.g. first element in array
# ${illuminarray[0]}
# ${minionarray[0]}

i=1
m=1

gunzip ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta.gz

pblat -tileSize=${TILE} -stepSize=${STEP} -threads=${NPROC} -noHead ${WORKDIR}${MINIONDIR}barcode${m}_all_pass.fasta ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl

gzip ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta

awk -v sim=${MATCH} '$1 >= sim { print }' ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl | gzip > ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl.gz

wc -l ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl.gz >> psl_wc-l.txt

# barcode 12 should have many good matches to A5 Papaver rhoeas
# barcode 12 should have 0 good matches to B5 Papaver somniferum (but is a congener of A5 so some false positives are expected)
# barcode 12 should have 0 good matches to A8 Tripleurospermum maritimum
# barcode 12 should have the most good matches to E3 Epilobium hirsutum

# the parallel version of the above.  creates lots of output files, which still need a command to be concatenated, and then need a way to select out only those Illumina reads that have > 220 matching nucleotides
# cat ${WORKDIR}${ILMNFAS}${QUERY}_bfc_trimed_uniques.fasta.gz | gunzip | ${PARALLEL} -j ${NPROC} --cat --pipe --block 10M --recstart ">" "blat -tileSize=${TILE} -stepSize=${STEP} -noHead ${WORKDIR}minION/barcode${TARGET}_all_pass.fasta {} ${WORKDIR}output_Plate1_${QUERY}_job{#}.psl" 2>${WORKDIR}output_Plate1_${QUERY}.stderr


# to create separate psl files for each minION read.  Not useful for my purposes
# mkdir output_Plate1_${QUERY}_output/
# cat output_Plate1_${QUERY}.psl | awk -v PFX=~/pollen_minION/output_Plate1_${QUERY}_output/ '{ file=PFX"/"$14".psl"; print $0>file; }'  # creates as many files as there are minION reference reads

ls blatout/
bpeek 207901


########## BSUB START ###################################################
#!/bin/sh
#BSUB -q short-eth  # short-eth has 24 hrs, mellanox-ib 7 days
#BSUB -n 24  # short-eth has 24 cores, add #BSUB -x # if i want an exclusive
#BSUB -J blatminse  # 10 chars max
#BSUB -B
#BSUB -N
#BSUB -oo blatminse.out
#BSUB -eo blatminse.err
#BSUB -R "rusage[mem=3800]" # 3800 = 3.8 GB RAM. BLAT uses < 500 MB
#BSUB -M 3800

. /etc/profile
# module load blat # v.36
PATH=$PATH:/gpfs/home/b042/scripts/pblat/ # pblat v.35  parallel blat
WORKDIR=~/pollen_minION/
MINIONDIR=minION/
ILMNFAS=plantskims/prepare_forBWA/prepare_forblast/
BLATDIR=blatout/
TILE=10
STEP=5
MATCH=213 # 85% of 251 is 213.  in other words, I am asking for the illumina read to have at least 85% similarity to the minION read
NPROC=$(grep -c processor /proc/cpuinfo)

cd ${WORKDIR}
mkdir ${BLATDIR}

declare -a illuminarray=(A1 A2 A3 B12 A4 A5 A12B5B6 A7A8 A9 A10 A11 B1 B2 C4 B3 B4 B7 B8 B9 B10 B11 C1 C2 C3 C5 C6C7 C8 C9 C10 C11 C12 D1 D2 D8 D3 D4 D5 D6 D7 D10 D9 D11 D12 E1 E2 E3 E4 E5 E6 E7 E8) # # Illumina reads = plant skims that have been adapter-trimmed, bfc-denoised, merged, dereplicated, and converted to fasta
declare -a minionarray=(01 02 03 04 06 08 09 10 11 12) # minION reads, all barcodes that successfully sequenced

rm -f psl_wc-l_Minsimilarity_${MATCH}.txt
touch psl_wc-l_Minsimilarity_${MATCH}.txt
echo "Illumina minION    Good_reads"  >> psl_wc-l_Minsimilarity_${MATCH}.txt

for i in "${illuminarray[@]}"
do
     for m in "${minionarray[@]}"
     do
          echo "Illumina reference ${i}"
          echo "minION query ${m}"

          gunzip ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta.gz

          pblat -tileSize=${TILE} -stepSize=${STEP} -threads=24 -noHead ${WORKDIR}${MINIONDIR}barcode${m}_all_pass.fasta ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl

          gzip ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta

          awk -v sim=${MATCH} '$1 >= sim { print }' ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl | gzip > ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl.gz

          rm ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl

          echo "${i}     ${m} $(wc -l ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl.gz)" >> psl_wc-l_Minsimilarity_${MATCH}.txt
     done
done
########## BSUB END ###################################################

# possibly could use anonymous named pipes to gzip psl output file

# (awk -v sim=${MATCH} '$1 >= sim { print }' | gzip > ) is applied to the output before being saved to disk

# blat -tileSize=${TILE} -stepSize=${STEP} -noHead ${WORKDIR}${MINIONDIR}barcode${m}_all_pass.fasta ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta.gz (awk -v sim=${MATCH} '$1 >= sim { print }' | gzip > ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl.gz)
#
# wc -l ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl.gz >> psl_wc-l.txt
