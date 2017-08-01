#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# shell script for BLAT mapping illumina reads to minION reads
#######################################################################################
#######################################################################################

# idea taken from NaS pipeline from Madoui et al. (2015, BMC Genomics)
# https://github.com/institut-de-genomique/NaS
# https://doi.org/10.1186/s12864-015-1519-z
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1519-z

# example commands taken from the NaS pipeline
              # blat -tileSize=$TILE -stepSize=$STEP -noHead barcode12 Plate1_A5_CCTCAGAGA_L002 output.psl
              # cat $OUTPUT_DIR/tmp/ILMN_reads.fa | parallel -j $NB_PROC --cat --pipe --block 10M --recstart ">" "$BLAT -tileSize=$TILE -stepSize=$STEP -noHead $NANO_READS {} $OUTPUT_DIR/tmp/psl/blat-alignment.job{#}.tile$TILE.step$STEP.psl" >$OUTPUT_DIR/tmp/blat-alignment.stderr
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

blat -tileSize=${TILE} -stepSize=${STEP} -noHead ${WORKDIR}minION/barcode12_all_pass.fasta ${WORKDIR}${ILMNFAS}A5_bfc_trimed_uniques.fasta.gz ${WORKDIR}output_Plate1_A5.psl

mkdir output_Plate1_A5_output/
cat output_Plate1_A5_cp.psl | awk -v PFX=~/pollen_minION/output_Plate1_A5_output/ '{ file=PFX"/"$14".psl"; print $0>file; }'  # creates as many files as there are minION reference reads


#### bsub version

#!/bin/sh
#BSUB -q mellanox-ib  # mellanox-ib 7 days
#BSUB -n 16  # add #BSUB -x # if i want an exclusive
#BSUB -J blatminion  # 10 chars max
#BSUB -B
#BSUB -N
#BSUB -oo blatminion.out
#BSUB -eo blatminion.err
#BSUB -R "rusage[mem=64000]" # 72000 = 64 GB RAM. mellanox-ib has 128 GB
#BSUB -M 73000


. /etc/profile
module load blat # v.36
PARALLEL=~/scripts/parallel-20170722/bin/parallel
ILMNFAS=plantskims/prepare_forBWA/prepare_forblast/
WORKDIR=~/pollen_minION/
TILE=10
STEP=5

blat -tileSize=${TILE} -stepSize=${STEP} -noHead ${WORKDIR}minION/barcode12_all_pass.fasta ${WORKDIR}${ILMNFAS}A5_bfc_trimed_uniques.fasta.gz ${WORKDIR}output_Plate1_A5.psl
