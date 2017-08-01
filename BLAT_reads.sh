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

# GNU parallel install on hpc
              # cd scripts/
              # wget http://ftp.heanet.ie/mirrors/gnu/parallel/parallel-latest.tar.bz2
              # bzip2 -d parallel-latest.tar.bz2
              # tar -xvf parallel-latest.tar
              # cd parallel-20170722
              # (wget -O - pi.dk/3 || curl pi.dk/3/ || fetch -o - http://pi.dk/3) | bash  # installs to ~/bin
              # mv ~/bin ~/scripts/parallel-20170722/
              # ~/scripts/parallel-20170722/bin/parallel -h
              # O. Tange (2011): GNU Parallel - The Command-Line Power Tool, ;login: The USENIX Magazine, February 2011:42-47.

# example commands taken from the NaS pipeline
              # blat -tileSize=$TILE -stepSize=$STEP -noHead barcode12 Plate1_A5_CCTCAGAGA_L002 output.psl
              # cat $OUTPUT_DIR/tmp/ILMN_reads.fa | parallel -j $NB_PROC --cat --pipe --block 10M --recstart ">" "$BLAT -tileSize=$TILE -stepSize=$STEP -noHead $NANO_READS {} $OUTPUT_DIR/tmp/psl/blat-alignment.job{#}.tile$TILE.step$STEP.psl" >$OUTPUT_DIR/tmp/blat-alignment.stderr
                            # This will BLAT map reads.  $TILE=10 and $STEP=5
              # mkdir output_dir
              # cat your_psl_file.psl | awk -v PFX=output_dir '{ file=PFX"/"$14".psl"; print $0>file; }'
                            # This will create one file per MinION read with the name of Illumina reads that mapped to this particular read.

#### Interactive version
ssh grace
ssh hpc.uea.ac.uk
interactive -x
# interactive -x -R "rusage[mem=20000]" -M 22000  # alternative is:  interactive -q interactive-lm
module load blat # v. 36
PARALLEL=~/scripts/parallel-20170722/bin/parallel
# ${PARALLEL} -h # to run
PARALLEL=~/scripts/parallel-20170722/bin/parallel
ILMNFAS=plantskims/prepare_forBWA/prepare_forblast/
WORKDIR=~/pollen_minION/
TILE=10
STEP=5

blat -tileSize=${TILE} -stepSize=${STEP} -noHead ${WORKDIR}minION/barcode12_all_pass.fasta ${WORKDIR}${ILMNFAS}A5_bfc_trimed_uniques.fasta.gz ${WORKDIR}output_Plate1_A5.psl


#### bsub version

#!/bin/sh
#BSUB -q huge-mem  # debug, short (5 hours), medium (24 hours), long (3 days)
#BSUB -n 16
#BSUB -x
#BSUB -J blatminion  # 10 chars max
#BSUB -B
#BSUB -N
#BSUB -oo blatminion.out
#BSUB -eo blatminion.err
#BSUB -R "rusage[mem=510000]"
#BSUB -M 512000


. /etc/profile
module load blat # v.36
PARALLEL=~/scripts/parallel-20170722/bin/parallel
ILMNFAS=plantskims/prepare_forBWA/prepare_forblast/
WORKDIR=~/pollen_minION/
TILE=10
STEP=5

blat -tileSize=${TILE} -stepSize=${STEP} -noHead ${WORKDIR}minION/barcode12_all_pass.fasta ${WORKDIR}${ILMNFAS}A5_bfc_trimed_uniques.fasta.gz ${WORKDIR}output_Plate1_A5.psl
