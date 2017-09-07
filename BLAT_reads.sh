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
# 1. parallel:  run all 50 skims against all 12 barcodes, try putting the loop into a shell script and calling that in the bsub file
# 2. script to separate psl files into one file for each minION read, and then filter by number of matching reads, evenness of reads, etc.?, and record in a table
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

# https://davetang.org/muse/2013/11/18/using-gnu-parallel/
#### Interactive version
ssh grace
ssh hpc.uea.ac.uk
# interactive -x # exclusive job, which is alternative to interactive -R
interactive -R "rusage[mem=4000]" -M 4000 # 4000=4GB RAM, which is the max for an hpc interactive job
module load blat # v. 36
# PATH=$PATH:/gpfs/home/b042/scripts/pblat/ # pblat v.35
PATH=$PATH:/gpfs/home/b042/scripts/parallel-20170722/bin/
WORKDIR=~/pollen_minION/
MINIONDIR=minION/
ILMNFAS=plantskims/prepare_forBWA/prepare_forblast/
BLATDIR=blatout_parallel_version/
TILE=10
STEP=5
MATCH=213 # 85% of 251 is 213.  in other words, I am asking for the illumina read to have at least 85% similarity to the minION read
# NPROC=$(grep -c processor /proc/cpuinfo)

cd ${WORKDIR}${BLATDIR}
# (A1 A2 A3 A4 A5 A8 A9 A11 B1 B2 B3     B4)
# (B5 B7 B8 B9 B10 B11 B12 C1 C2 C3 C4 C5)
# (C6 C8 C10 C11 C12 D1 D2 D8 D3 D4 D5 D6)
# (D7 D9 D10 D11 D12 E1 E2 E3 E4 E5 E6 E7 E8)

cd ${WORKDIR}${ILMNFAS}
# (A1 A2 A3 A4 A5 A8 A9 A11 B1 B2 B3 B4 B5 B7 B8 B9 B10 B11 B12 C1 C2 C3 C4 C5 C6 C8 C10 C11 C12 D1 D2 D8 D3 D4 D5 D6 D7 D9 D10 D11 D12 E1 E2 E3 E4 E5 E6 E7 E8)
# (01 02 03 04 06 08 09 10 11 12)
i=A5
m=01

# gunzip ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta.gz # blat requires fasta file input

# Biostars code, which automatically cats the partial outputs before redirecting output to a single *.psl file
# cat foo.fa | parallel --round-robin --pipe --recstart '>' 'blat -noHead genome.fa stdin >(cat) >&2' >foo.psl
# this syntax is from the NaS's own code.
# cat $OUTPUT_DIR/tmp/ILMN_reads.fa | parallel -j $NB_PROC --block 10M --pipe --recstart ">" "$BLAT -tileSize=$TILE -stepSize=$STEP -noHead $NANO_READS stdin >(cat) >&2" > $OUTPUT_DIR/tmp/blat-alignment.tile$TILE.step$STEP.psl 2>$OUTPUT_DIR/tmp/blat-alignment.stderr

# My version of the above.
# try gunzipping and piping fasta.gz file. then shouldn't need to gzip afterwards
# this works
nohup zcat ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta.gz | parallel -j 16 --round-robin --cat --pipe --block 10M --recstart ">" "blat -tileSize=${TILE} -stepSize=${STEP} -noHead ${WORKDIR}${MINIONDIR}barcode${m}_all_pass.fasta {} >(cat) >&2" > ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl 2> ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.stderr &


# NaS author-suggested code to me. This works but needs to be catted afterwards
# nohup gunzip ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta | cat | parallel -j 16 --cat --pipe --block 10M --recstart ">" "blat -tileSize=${TILE} -stepSize=${STEP} -noHead ${WORKDIR}${MINIONDIR}barcode${m}_all_pass.fasta {} ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}_job{#}.psl" > ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}_sepjobs.stderr &
# blat {options} minionreads illuminareads output.psl
     # -j 16 means 16 jobs at a time = number of cores available
     # --cat make a temp file with filename in {}, so {} can be used as input
     # --pipe pipe contents of stdin instead of the using the stdin as arguments. --cat and --pipe are used together
     # --block 10M read 10M of data at a time
     # --recstart ">" each record starts with a ">", b/c is a fasta file. thus, won't cut a record in half
     # --round-robin only create the number of commands as number of jobs
     # {#} is the job number
     # {} is the stdin block from the illumina reads file

# cat ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}_job*.psl > ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl

# rm -f ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}_job*.psl
# rm -f output_ILMNREF_A11_MINIONBARCODE_01_sepjobs.stderr

# This will create one file per MinION read with the name of Illumina reads that mapped to this particular read. $14 is the column with the name of the minION read. $0 means the whole line.
# cat $OUTPUT_DIR/tmp/blat-alignment.tile$TILE.step$STEP.psl | awk -v PFX=$OUTPUT_DIR/reads '{ file=PFX"/"$14".psl"; print $0>file; }'
# cat your_psl_file.psl | awk -v PFX=output_dir '{ file=PFX"/"$14".psl"; print $0>file; }'

gzip ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta

awk -v sim=${MATCH} '$1 >= sim { print }' ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl | gzip > ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}_min_${MATCH}.psl.gz

# rm -f ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl

echo -e "${i}\t${m}\t$(wc -l ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl.gz)" >> psl_wc-l_Minsimilarity_${MATCH}.txt # tab delimited, -e allows echo to understand \t

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
bpeek 208877
head psl_wc-l_Minsimilarity_213.txt

# version with pblat
########## BSUB START ###################################################
#!/bin/sh
#BSUB -q mellanox-ib  # short-eth has 24 hrs, mellanox-ib 7 days
#BSUB -n 24  # mellanox-ib has 24 cores, add #BSUB -x # if i want an exclusive
#BSUB -J blatA1E8  # 10 chars max
#BSUB -B
#BSUB -N
#BSUB -oo blatA1E8.out
#BSUB -eo blatA1E8.err
#BSUB -R "rusage[mem=3800]" # 3800 = 3.8 GB RAM. BLAT uses < 500 MB
#BSUB -M 3800

. /etc/profile
module load blat # v.36
PATH=$PATH:/gpfs/home/b042/scripts/pblat/ # pblat v.35  parallel blat
WORKDIR=~/pollen_minION/
MINIONDIR=minION/
ILMNFAS=plantskims/prepare_forBWA/prepare_forblast/
BLATDIR=blatout/
TILE=10
STEP=5
MATCH=213 # 85% of 251 is 213.  in other words, I am asking for the illumina read to have at least 85% similarity to the minION read
NPROC=24 # 24 for mellanox-ib, 24 for long-eth

cd ${WORKDIR}  # if cd fails
# if starting new, uncomment the following four lines
# mkdir ${BLATDIR}
# rm -f psl_wc-l_Minsimilarity_${MATCH}.txt
touch ${WORKDIR}${BLATDIR}psl_wc-l_Minsimilarity_${MATCH}.txt
echo -e "Illumina\tminION\tGood_reads"  >> ${WORKDIR}${BLATDIR}psl_wc-l_Minsimilarity_${MATCH}.txt

# A7, A12, B6 and C7 are duplicates of other skims
# A6, A10 and C9 did not produce reference libraries (A10 produced a very small one)
declare -a illuminarray=(A1 A2 A3 A4 A5 A8 A9 A11 B1 B2 B3 B4 B5 B7 B8 B9 B10 B11 B12 C1 C2 C3 C4 C5 C6 C8 C10 C11 C12 D1 D2 D8 D3 D4 D5 D6 D7 D9 D10 D11 D12 E1 E2 E3 E4 E5 E6 E7 E8) # Illumina reads = plant skims that have been adapter-trimmed, bfc-denoised, merged, dereplicated, and converted to fasta
declare -a minionarray=(01 02 03 04 06 08 09 10 11 12) # minION reads, all barcodes that successfully sequenced

# smaller test set
# declare -a illuminarray=(A2 A9) # # Illumina reads = plant skims that have been adapter-trimmed, bfc-denoised, merged, dereplicated, and converted to fasta
# declare -a minionarray=(01 02) # minION reads, all barcodes that successfully sequenced


for i in "${illuminarray[@]}"
do
     gunzip ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta.gz

     for m in "${minionarray[@]}"
     do
          echo "Illumina reference ${i}"
          echo "minION query ${m}"

          pblat -tileSize=${TILE} -stepSize=${STEP} -threads=${NPROC} -noHead ${WORKDIR}${MINIONDIR}barcode${m}_all_pass.fasta ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl

          # this version works when input singly but within the loop
          # zcat ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta.gz | parallel -j ${NPROC} --round-robin --cat --pipe --block 10M --recstart ">" "blat -tileSize=${TILE} -stepSize=${STEP} -noHead ${WORKDIR}${MINIONDIR}barcode${m}_all_pass.fasta {} >(cat) >&2" > ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl 2> ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.stderr
          # blat {options} minionreads illuminareads output.psl
               # -j ${NPROC} means NPROC jobs at a time = number of cores available
               # --cat make a temp file with filename in {}, so {} can be used as input
               # --pipe pipe contents of stdin instead of the using the stdin as arguments. --cat and --pipe are used together
               # --block 10M read 10M of data at a time
               # --recstart ">" each record starts with a ">", b/c is a fasta file. thus, won't cut a record in half
               # --round-robin only create the number of commands as number of jobs
               # {} is the stdin block from the illumina reads file

          awk -v sim=${MATCH} '$1 >= sim { print }' ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl > ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}_min_${MATCH}.psl

          rm -f ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl

          echo -e "${i}\t${m}\t$(wc -l ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}_min_${MATCH}.psl)" >> ${WORKDIR}${BLATDIR}psl_wc-l_Minsimilarity_${MATCH}.txt # tab delimited, -e allows echo to understand \t
     done

     gzip ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta

done
########## BSUB END ###################################################

# was not able to use anonymous named pipes to gunzip the input fasta.gz file
# possibly could use anonymous named pipes to awk and gzip psl output file

# pblat -tileSize=${TILE} -stepSize=${STEP} -threads=24 -noHead ${WORKDIR}${MINIONDIR}barcode${m}_all_pass.fasta ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta >(awk -v sim=${MATCH} '$1 >= sim { print }' ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl | gzip > ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl.gz) is applied to the output before being saved to disk

# create multiple files from a single file plus an index in the filename
tee blatminion_short-eth_20170826_{01,02,03,04,06,08,09,10,11,12}.bsub < blatminion_short-eth_20170825.bsub >/dev/null

# to submit multiple bsub files
cd ~/pollen_minION
declare -a minionarray=(01 02 03 04 06 08 09 10 11 12)
for m in "${minionarray[@]}"
do
     bsub < blatminion_short-eth_20170826_${m}.bsub
done

bjobs
#
ls ${WORKDIR}${ILMNFAS}*_bfc_trimed_uniques.fasta # to see which Illumina query files are being blatted
ls ${WORKDIR}${ILMNFAS}*_bfc_trimed_uniques.fasta.gz # to see gzipped Illumina query files
ls -l blatout/output_ILMNREF_A4_MINIONBARCODE*.psl.gz


# blatall.bsub, zcat to parallel to cat
# run in pollen_minION/
########## BSUB START ###################################################
#!/bin/sh
#BSUB -q mellanox-ib  # short-ib has 24 hrs, mellanox-ib 7 days
#BSUB -n 28  # mellanox-ib has 28 cores, add #BSUB -x # if i want an exclusive
#BSUB -J blatall  # 10 chars max
#BSUB -B
#BSUB -N
#BSUB -oo blatall.out
#BSUB -eo blatall.err
#BSUB -R "rusage[mem=12000]" # 12000 = 12 GB RAM. BLAT typically uses < 5000 MB
#BSUB -M 12000

. /etc/profile
module load blat # v.36

bash blatall.sh
########## BSUB END ###################################################


########## blatall.sh, place in pollen_minION/ ###################################################
#!/bin/bash
set -e
set -u
set -o pipefail

WORKDIR=~/pollen_minION/
MINIONDIR=minION/
ILMNFAS=plantskims/prepare_forBWA/prepare_forblast/
BLATDIR=blatout_parallel_version/
TILE=10
STEP=5
MATCH=213 # 85% of 251 is 213.  in other words, I am asking for the illumina read to have at least 85% similarity to the minION read
NPROC=28 # 28 for mellanox-ib, 24 for long-eth

cd ${WORKDIR}
# if starting new, uncomment the following four lines
# mkdir ${BLATDIR}
# rm -f psl_wc-l_Minsimilarity_${MATCH}.txt
# touch ${WORKDIR}${BLATDIR}psl_wc-l_Minsimilarity_${MATCH}.txt
# echo -e "Illumina\tminION\tGood_reads"  >> ${WORKDIR}${BLATDIR}psl_wc-l_Minsimilarity_${MATCH}.txt

# A7, A12, B6 and C7 are duplicates of other skims
# A6, A10 and C9 did not produce reference libraries (A10 produced a very small one)
# in practice, i split the whole array into four and ran separate bsub jobs.
# see below for how i processed the huge A5-01 combination
declare -a illuminarray=(A1 A2 A3 A4 A5 A8 A9 A11 B1 B2 B3 B4 B5 B7 B8 B9 B10 B11 B12 C1 C2 C3 C4 C5 C6 C8 C10 C11 C12 D1 D2 D3 D4 D5 D6 D7 D8 D9 D10 D11 D12 E1 E2 E3 E4 E5 E6 E7 E8) # Illumina reads = 49 plant skims that have been adapter-trimmed, bfc-denoised, merged, dereplicated, and converted to fasta
declare -a minionarray=(01 02 03 04 06 08 09 10 11 12) # 10 minION read files, all barcodes that successfully sequenced

# test set with small files
# declare -a illuminarray=(A2 B11) # # Illumina reads = plant skims that have been adapter-trimmed, bfc-denoised, merged, dereplicated, and converted to fasta
# declare -a minionarray=(02 06) # minION reads, all barcodes that successfully sequenced

for i in "${illuminarray[@]}"
do
     for m in "${minionarray[@]}"
     do
          echo "Starting blat for Illumina reference ${i} and minION query ${m}.  $(date)"

          zcat ${WORKDIR}${ILMNFAS}${i}_bfc_trimed_uniques.fasta.gz | parallel -j ${NPROC} --round-robin --cat --pipe --block 10M --recstart ">" "blat -tileSize=${TILE} -stepSize=${STEP} -noHead ${WORKDIR}${MINIONDIR}barcode${m}_all_pass.fasta {} >(cat) >&2" > ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl 2> ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.stderr
          # blat {options} minionreads illuminareads output.psl
               # -j ${NPROC} means NPROC jobs at a time = number of cores available
               # --cat make a temp file with filename in {}, so {} can be used as input
               # --pipe pipe contents of stdin instead of the using the stdin as arguments. --cat and --pipe are used together
               # --block 10M read 10M of data at a time
               # --recstart ">" each record starts with a ">", b/c input is a fasta file. thus, won't cut a record in half
               # --round-robin only create the number of commands as number of jobs
               # {} is the stdin block from the illumina reads file

          echo "Finished blat for Illumina reference ${i} and minION query ${m}.  $(date)"

          awk -v sim=${MATCH} '$1 >= sim { print }' ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl > ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}_min_${MATCH}.psl

          rm -f ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}.psl

          echo -e "${i}\t${m}\t$(wc -l ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}_min_${MATCH}.psl)" >> ${WORKDIR}${BLATDIR}psl_wc-l_Minsimilarity_${MATCH}.txt # tab delimited, -e allows echo to understand \t

          echo -e "${i}\t${m}\t$(wc -l ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}_min_${MATCH}.psl)"
     done
done
########## blatall.sh end ###################################################

# For the A5-01 combination, the above parallel version of the command crashed (huge files created), so i used this version, which creates multiple output files and then catted them
     # zcat ${WORKDIR}${ILMNFAS}A5_bfc_trimed_uniques.fasta.gz | parallel -j ${NPROC} --round-robin --cat --pipe --block 10M --recstart ">" "blat -tileSize=${TILE} -stepSize=${STEP} -noHead ${WORKDIR}${MINIONDIR}barcode01_all_pass.fasta {} ${WORKDIR}${BLATDIR}output_ILMNREF_A5_MINIONBARCODE_01_job{#}.psl"
     # however, this meant that i had to awk each job output separately (which i also did with parallel) and then i catted the outputs
          # "ls" -1 output_ILMNREF_A5_MINIONBARCODE_01_job*.psl > A501filestoawk # create list of files to awk. quoting "ls" overrides the alias that i set up
          # cat A501filestoawk # confirm list
          # my_awk='$1 >= sim { print }' # store awk command in a variable because it can conflict with parallel
          #
          # # output filename is input filename minus the extension, if i use {.}
          # # in interactive job mode, this runs 16 parallel jobs
          # # -dryrun to see commands before i run them, -k to keep filenames in order (helpful for debugging)
          # parallel --dryrun -k "cat {} | awk -v sim=${MATCH} '$my_awk' > {.}_min_${MATCH}.psl" :::: A501filestoawk
          # # live command
          # parallel "cat {} | awk -v sim=${MATCH} '$my_awk' > {.}_min_${MATCH}.psl" :::: A501filestoawk
          #
          # cat output_ILMNREF_A5_MINIONBARCODE_01_job{1..24}_min_213.psl > output_ILMNREF_A5_MINIONBARCODE_01_min_213.psl
          #
          # rm output_ILMNREF_A5_MINIONBARCODE_01_job{1..24}.psl
          # rm output_ILMNREF_A5_MINIONBARCODE_01_job{1..24}_min_213.psl
          # rm A501filestoawk
          # echo -e "A5\t01\t$(wc -l ${WORKDIR}${BLATDIR}output_ILMNREF_A5_MINIONBARCODE_01_min_213.psl)" >> ${WORKDIR}${BLATDIR}psl_wc-l_Minsimilarity_${MATCH}.txt
          # gzip output_ILMNREF_A5_MINIONBARCODE_01_min_213.psl



########## interactive awkblat for different match values ###################################################
WORKDIR=~/pollen_minION/
BLATDIR=blatout_parallel_version/

MATCH=225 # 90% of 251 is 225.5. keep only illumina reads with at least 90% similarity to the minION read
cd ${WORKDIR}${BLATDIR}
#
# A7, A12, B6 and C7 are duplicates of other skims
# A6, A10 and C9 did not produce reference libraries (A10 produced a very small one)
declare -a illuminarray=(A1 A2 A3 A4 A5 A8 A9 A11 B1 B2 B3 B4 B5 B7 B8 B9 B10 B11 B12 C1 C2 C3 C4 C5 C6 C8 C10 C11 C12 D1 D2 D3 D4 D5 D6 D7 D8 D9 D10 D11 D12 E1 E2 E3 E4 E5 E6 E7 E8) # Illumina reads = 49 plant skims that have been adapter-trimmed, bfc-denoised, merged, dereplicated, and converted to fasta
declare -a minionarray=(01 02 03 04 06 08 09 10 11 12) # 10 minION read files, all barcodes that successfully sequenced

my_awk='$1 >= sim { print }'
parallel "echo "Starting gunzip for Illumina reference {1} and minION query {2}."; gunzip ${WORKDIR}${BLATDIR}output_ILMNREF_{1}_MINIONBARCODE_{2}_min_213.psl.gz; echo "Starting awk for Illumina reference {1} and minION query {2}."; awk -v sim=${MATCH} '${my_awk}' ${WORKDIR}${BLATDIR}output_ILMNREF_{1}_MINIONBARCODE_{2}_min_213.psl > ${WORKDIR}${BLATDIR}output_ILMNREF_{1}_MINIONBARCODE_{2}_min_${MATCH}.psl; echo "Starting gzip for Illumina reference {1} and minION query {2}."; gzip ${WORKDIR}${BLATDIR}output_ILMNREF_{1}_MINIONBARCODE_{2}_min_213.psl" ::: "${illuminarray[@]}" ::: "${minionarray[@]}"

# check
ls ${WORKDIR}${BLATDIR}output_ILMNREF_*_MINIONBARCODE_*_min_${MATCH}.psl | wc -l  # should be 490

# original loop version
# for i in "${illuminarray[@]}"
# do
#      for m in "${minionarray[@]}"
#      do
#           echo "Starting gunzip for Illumina reference ${i} and minION query ${m}.  $(date)"
#           gunzip ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}_min_213.psl.gz
#
#           echo "Starting awk for Illumina reference ${i} and minION query ${m}.  $(date)"
#           awk -v sim=${MATCH} '$1 >= sim { print }' ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}_min_213.psl > ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}_min_${MATCH}.psl
#
#           echo "Starting gzip for Illumina reference ${i} and minION query ${m}.  $(date)"
#           gzip ${WORKDIR}${BLATDIR}output_ILMNREF_${i}_MINIONBARCODE_${m}_min_213.psl
#      done
# done


########## blatall.sh end ###################################################


# interactive bedops and bedtools
ssh grace
ssh hpc.uea.ac.uk
interactive -R "rusage[mem=4000]" -M 4000 # 4000=4GB RAM, which is the max for an hpc interactive job
module load blat # v. 36
module load samtools # v. 1.5  for samtools faidx
module load java
module load bcftools
module load python
PATH=$PATH:/gpfs/home/b042/scripts/parallel-20170722/bin/
PATH=$PATH:/gpfs/home/b042/scripts/bedbin/
PATH=$PATH:/gpfs/home/b042/scripts/bedtools2/bin/
PATH=$PATH:/gpfs/home/b042/scripts/depth-cover/depth-cover.jar
WORKDIR=~/pollen_minION/
MINIONDIR=minION/
ILMNFAS=plantskims/prepare_forBWA/prepare_forblast/
BLATDIR=blatout_parallel_version/
TILE=10
STEP=5
MATCH=225 # 90% of 251 is 225.  in other words, I am asking for the illumina read to have at least 85% similarity to the minION read

cd ${WORKDIR}${BLATDIR}

# # just do once:  create "genome" file of the minION reads, which are the "chromosomes" expected by bedtools
# # https://www.biostars.org/p/70795/
# # generate summary faidx files
# cd ${WORKDIR}${MINIONDIR}
# "ls" -1 barcode*_all_pass.fasta > filestofaidx
# cat filestofaidx
# parallel --dryrun "samtools faidx {}" :::: filestofaidx
# parallel "samtools faidx {}" :::: filestofaidx
# rm filestofaidx
# # keep only first two columns of the faidx files
# "ls" -1 barcode*_all_pass.fasta.fai > filestoawk
# cat filestoawk
# my_awk='{ print $1,$2 }'
# parallel --dryrun "awk -v OFS='\t' '${my_awk}' {} > {.}.txt" :::: filestoawk
# parallel "awk -v OFS='\t' '${my_awk}' {} > {.}.txt" :::: filestoawk
# rm -f filestoawk

# convert psl to bed
cd ${WORKDIR}${BLATDIR}
"ls" -1 output_ILMNREF_*_MINIONBARCODE_*${MATCH}.psl > filestobed
cat filestobed # | wc -l # should be 490 files
parallel "cat {} | psl2bed - | gzip > {.}.bed.gz" :::: filestobed
rm filestobed

# convert bed to bam
cd ${WORKDIR}${BLATDIR}
declare -a illuminarray=(A1 A2 A3 A4 A5 A8 A9 A11 B1 B2 B3 B4 B5 B7 B8 B9 B10 B11 B12 C1 C2 C3 C4 C5 C6 C8 C10 C11 C12 D1 D2 D3 D4 D5 D6 D7 D8 D9 D10 D11 D12 E1 E2 E3 E4 E5 E6 E7 E8) # Illumina reads = 49 plant skims that have been adapter-trimmed, bfc-denoised, merged, dereplicated, and converted to fasta
declare -a minionarray=(01 02 03 04 06 08 09 10 11 12) # 10 minION read files, all barcodes that successfully sequenced
# bedToBam -i rmsk.hg18.chr21.bed -g human.hg18.genome > rmsk.hg18.chr21.bam
# changed bed to bam format and sort output
parallel "bedToBam -i ${WORKDIR}${BLATDIR}output_ILMNREF_{1}_MINIONBARCODE_{2}_min_${MATCH}.bed.gz -g ${WORKDIR}${MINIONDIR}barcode{2}_all_pass.fasta.txt > ${WORKDIR}${BLATDIR}output_ILMNREF_{1}_MINIONBARCODE_{2}_min_${MATCH}.bam | samtools sort -o ${WORKDIR}${BLATDIR}output_ILMNREF_{1}_MINIONBARCODE_{2}_min_${MATCH}.bam" ::: "${illuminarray[@]}" ::: "${minionarray[@]}"
# check
ls ${WORKDIR}${BLATDIR}output_ILMNREF_*_MINIONBARCODE_*_min_${MATCH}.bam | wc -l  # should be 490



## calculate genome coverage bedtools

# for now, do this for subsets in test/ folder
# mkdir test
ILMN="A3"
cd ${WORKDIR}${BLATDIR}
# cp ${WORKDIR}${BLATDIR}output_ILMNREF_${ILMN}_MINIONBARCODE_*_min_${MATCH}.bam ${WORKDIR}${BLATDIR}test/
cp ${WORKDIR}${BLATDIR}output_ILMNREF_${ILMN}_MINIONBARCODE_*_min_${MATCH}.bed.gz ${WORKDIR}${BLATDIR}test/
cd test
ls

declare -a illuminarray=(${ILMN}) # Illumina reads = 49 plant skims that have been adapter-trimmed, bfc-denoised, merged, dereplicated, and converted to fasta
declare -a minionarray=(01 02 03 04 06 08 09 10 11 12) # 10 minION read files, all barcodes that successfully sequenced

# bedtools genomecov -d -i exons.bed -g genome.txt
# -d  # calculates coverage for all positions
parallel "bedtools genomecov -d -i ${WORKDIR}${BLATDIR}test/output_ILMNREF_{1}_MINIONBARCODE_{2}_min_${MATCH}.bed.gz -g ${WORKDIR}${MINIONDIR}barcode{2}_all_pass.fasta.txt | gzip > output_ILMNREF_{1}_MINIONBARCODE_{2}_min_${MATCH}.bed.d.genomecov.gz" ::: "${illuminarray[@]}" ::: "${minionarray[@]}"

zcat output_ILMNREF_${ILMN}_MINIONBARCODE_01_min_225.bed.d.genomecov.gz | less

# bg
parallel "bedtools genomecov -i ${WORKDIR}${BLATDIR}test/output_ILMNREF_{1}_MINIONBARCODE_{2}_min_${MATCH}.bed.gz -g ${WORKDIR}${MINIONDIR}barcode{2}_all_pass.fasta.txt > output_ILMNREF_{1}_MINIONBARCODE_{2}_min_${MATCH}.bed.bg.genomecov" ::: "${illuminarray[@]}" ::: "${minionarray[@]}"

cat output_ILMNREF_${ILMN}_MINIONBARCODE_01_min_225.bed.bg.genomecov | less


rm -f reads_bga.txt
rm -f reads_bg.txt
for i in "${minionarray[@]}"
do
     cut -f 1 output_ILMNREF_B2_MINIONBARCODE_${i}_min_${MATCH}.bed.bg.genomecov | uniq | wc -l >> reads_bg.txt
     cut -f 1 output_ILMNREF_B2_MINIONBARCODE_${i}_min_${MATCH}.bed.bga.genomecov | uniq | wc -l >> reads_bga.txt
done
cat reads_bga.txt
cat reads_bg.txt

### START HERE:  need to figure out a way to parse genome coverage and evenness from the d.genomecov files





#### practice code.
# samtools view -H output_ILMNREF_B2_MINIONBARCODE_12_min_225.bam | less
# grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'q
#
# samtools depth -a output_ILMNREF_B2_MINIONBARCODE_12_min_225.sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}'
#
# parallel "samtools depth output_ILMNREF_{1}_MINIONBARCODE_{2}_min_${MATCH}.bam | awk  '{sum+=$3} END { print "Average = ",sum/NR}' > " ::: "${illuminarray[@]}" ::: "${minionarray[@]}"
#
# samtools depth  *bamfile*  |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
#
# samtools depth *bamfile*  |  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}'












less output_ILMNREF_B2_MINIONBARCODE_01_min_238.bed.bga.genomecov

# Plan for R script
https://www.biostars.org/p/5165/  # samtools for calculating mean coverage and stdev (or evenness)
# convert bed to bam
bedtools bedtobam -bed12 -i input.bed -g hg38.chrom.sizes > output.bam


# for each psl file, convert to bed file:
#      http://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/psl2bed.html#downloads
#      http://bedtools.readthedocs.io/en/latest/#
#      https://github.com/bsipos/uncle_psl
#    calculate number of sequences of the illumina file, for normalisation
#    run an R script to do the following
#      group_by minION read, create a table for each minion read with
#      number of illumina reads
#      length of minion read
#      evenness of illumina read start positions
#      calculate number of sequences of the minION file, for normalisation
#      add number of sequences from the illumina file
#      calculate an evenness score for the illumina reads on that minion read
#
#    done
# done
#
