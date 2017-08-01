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

# example commands taken from the NaS pipeline
              # cat $OUTPUT_DIR/tmp/ILMN_reads.fa | parallel -j $NB_PROC --cat --pipe --block 10M --recstart ">" "$BLAT -tileSize=$TILE -stepSize=$STEP -noHead $NANO_READS {} $OUTPUT_DIR/tmp/psl/blat-alignment.job{#}.tile$TILE.step$STEP.psl" >$OUTPUT_DIR/tmp/blat-alignment.stderr
                            # This will BLAT map reads.  $TILE=10 and $STEP=5
              # mkdir output_dir
              # cat your_psl_file.psl | awk -v PFX=output_dir '{ file=PFX"/"$14".psl"; print $0>file; }'
                            # This will create one file per MinION read with the name of Illumina reads that mapped to this particular read.


ssh grace
ssh hpc.uea.ac.uk
interactive
module load blat # v. 36
PARALLEL=~/scripts/parallel-20170722/bin/parallel
# ${PARALLEL} -h # to run
