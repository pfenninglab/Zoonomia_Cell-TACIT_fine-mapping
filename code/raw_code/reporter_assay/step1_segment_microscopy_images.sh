#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen3,gpu
#SBATCH --time=24:00:00
#SBATCH --job-name=image
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=10G
#SBATCH --error=logs/segment_nuclei_%A_%a.txt
#SBATCH --output=logs/segment_nuclei_%A_%a.txt
#SBATCH --array=1-112%4

source activate cp4

PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
DATADIR=${PROJDIR}/data/raw_data/reporter_assay
CODEDIR=${PROJDIR}/code/raw_code/reporter_assay

cd $CODEDIR

FILE_NAME=$(ls ${DATADIR}/images/*/*.nd2 | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo $(basename ${FILE_NAME}) 
python ${CODEDIR}/process_imgs.py --file "${FILE_NAME}" \
--prefix 'Batch1' --out-dir ${DATADIR}
