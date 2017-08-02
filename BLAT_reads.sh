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
module load blat # v. 36
PARALLEL=~/scripts/parallel-20170722/bin/parallel
ILMNFAS=plantskims/prepare_forBWA/prepare_forblast/
WORKDIR=~/pollen_minION/
TILE=10
STEP=5
TARGET=12  # minION reads
QUERY=E3  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta
NPROC=$(grep -c processor /proc/cpuinfo)

cd ${WORKDIR}; blat -tileSize=${TILE} -stepSize=${STEP} -noHead ${WORKDIR}minION/barcode${TARGET}_all_pass.fasta ${WORKDIR}${ILMNFAS}${QUERY}_bfc_trimed_uniques.fasta.gz ${WORKDIR}output_Plate1_${QUERY}.psl
# barcode 12 should have many good matches to A5 Papaver rhoeas
# barcode 12 should have 0 good matches to B5 Papaver somniferum (but is a congener of A5 so some false positives are expected)
# barcode 12 should have 0 good matches to A8 Tripleurospermum maritimum
# barcode 12 should have the most good matches to E3 Epilobium hirsutum

# the parallel version of the above.  creates lots of output files, which still need a command to be concatenated, and then need a way to select out only those Illumina reads that have > 220 matching nucleotides
# cat ${WORKDIR}${ILMNFAS}${QUERY}_bfc_trimed_uniques.fasta.gz | gunzip | ${PARALLEL} -j ${NPROC} --cat --pipe --block 10M --recstart ">" "blat -tileSize=${TILE} -stepSize=${STEP} -noHead ${WORKDIR}minION/barcode${TARGET}_all_pass.fasta {} ${WORKDIR}output_Plate1_${QUERY}_job{#}.psl" 2>${WORKDIR}output_Plate1_${QUERY}.stderr


# to create separate psl files for each minION read.  Not useful for my purposes
# mkdir output_Plate1_${QUERY}_output/
# cat output_Plate1_${QUERY}.psl | awk -v PFX=~/pollen_minION/output_Plate1_${QUERY}_output/ '{ file=PFX"/"$14".psl"; print $0>file; }'  # creates as many files as there are minION reference reads


#### bsub version

#!/bin/sh
#BSUB -q short-eth  # short-eth has 24 hrs, mellanox-ib 7 days
#BSUB -n 12  # short-eth has 24 cores, add #BSUB -x # if i want an exclusive
#BSUB -J blatmin  # 10 chars max
#BSUB -B
#BSUB -N
#BSUB -oo blatmin.out
#BSUB -eo blatmin.err
#BSUB -R "rusage[mem=3000]" # 3000 = 3 GB RAM. BLAT uses < 500 MB
#BSUB -M 4000


. /etc/profile
module load blat # v.36
PARALLEL=~/scripts/parallel-20170722/bin/parallel
ILMNFAS=plantskims/prepare_forBWA/prepare_forblast/
WORKDIR=~/pollen_minION/
TILE=10
STEP=5
TARGET=12  # minION reads
QUERY=E3  # Illumina reads = plant skims that have been adapter-trimmed, denoised, merged, and converted to fasta

blat -tileSize=${TILE} -stepSize=${STEP} -noHead ${WORKDIR}minION/barcode${TARGET}_all_pass.fasta ${WORKDIR}${ILMNFAS}${QUERY}_bfc_trimed_uniques.fasta.gz ${WORKDIR}output_Plate1_${QUERY}.psl
