#!/bin/bash
#SBATCH --partition=pool1
#SBATCH --time 1-0
#SBATCH --job-name=scOpen
#SBATCH --mem=40G
#SBATCH --error=logs/scOpen_%A_%a.out.txt
#SBATCH --output=logs/scOpen_%A_%a.out.txt
#SBATCH --array 1

# state the log file.
PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/raw_data/cross_species_peak_orthologs
cd $PROJDIR

# find the sparse matrix file
INPUT=`ls -d ${PROJDIR}/scOpen/*_OrthologPeakMatrix| awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = "," } ;NR==IND {print $1}'`
OUTPUT_PREFIX=$(basename ${INPUT})

echo 'Performing scOpen on '${OUTPUT_PREFIX}'.'

# run scOpen 
scopen --input ${INPUT} \
--output-dir ${PROJDIR}/scOpen \
--output-prefix ${OUTPUT_PREFIX} \
--input-format 10X --output-format dense \
--n-components 30 --max-iter 1000 --verbose 1
