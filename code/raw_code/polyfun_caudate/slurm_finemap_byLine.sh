#!/bin/bash
#SBATCH --partition=pfen1
#SBATCH --mem=60G
#SBATCH --job-name=runLine
#SBATCH --error=logs/runLine_%A_%a.txt
#SBATCH --output=logs/runLine_%A_%a.txt

source ~/.bashrc; conda activate polyfun;

# just in case array wildcard '@' is in the JobsFileName
JobsFileName=$(echo $JobsFileName | sed "s/@/$SLURM_ARRAY_TASK_ID/g")
NumJobs=$(wc -l $JobsFileName | cut -d ' ' -f1)
export TEMP=/scratch/bnphan/polyfun; mkdir -p $TEMP

# switch out LD matrix. don't download to tmp
CACHEDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/polyfun/LD_cache
CACHEDIR1=$(echo $CACHEDIR | sed 's_/_\\/_g')
LDDIR=https://data.broadinstitute.org/alkesgroup/UKBB_LD
LDDIR1=$(echo $LDDIR | sed 's_/_\\/_g')

for ID in `seq $NumJobs`; do 
echo -e "Reading in $(basename ${JobsFileName}) and running line ${ID} of ${NumJobs}."
# grab the line to run
RUNME=$(cat ${JobsFileName} | awk "NR==${ID}")

# make sure LD block has been downloaded to cachedir 
LDFILE=$(echo $RUNME | awk -F"--ld" '{print $2}'| awk -F" " '{print $1}')
if [ ! -f $CACHEDIR/$(basename $LDFILE).gz ]; then wget ${LDFILE}.gz -P $CACHEDIR; fi
if [ ! -f $CACHEDIR/$(basename $LDFILE).npz ]; then wget ${LDFILE}.npz -P $CACHEDIR; fi

# change LD matrix location to cache dir
RUNME=$(echo $RUNME | sed "s/$LDDIR1/$CACHEDIR1/g" )

# extract the output file
OUTFILE=$(echo $RUNME | awk -F"--out" '{print $2}'| awk -F" " '{print $1}')
gzip -t $OUTFILE; # remove if gzipped file is corrupt
if [[ $(echo $?) != '0' ]]; then rm $OUTFILE; fi
if [[ ! -f $OUTFILE ]]; then echo "Finemap region not found, generating: $(basename ${OUTFILE})"; $RUNME; fi
done

