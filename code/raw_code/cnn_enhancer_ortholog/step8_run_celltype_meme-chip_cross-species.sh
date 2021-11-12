#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --time 3-0
#SBATCH --job-name=meme-chip
#SBATCH --mem=10G
#SBATCH --error=logs/meme-chip_%A_%a_out.txt
#SBATCH --output=logs/meme-chip_%A_%a_out.txt
#SBATCH --array=1-24

SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/cnn_enhancer_ortholog
DATADIR=${SETDIR}/data/raw_data/cnn_enhancer_ortholog
cd $CODEDIR

#############################################
# get the cell type to be used for training #
CELLS=( NULL MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Microglia Oligo OPC )
(( INDC=($SLURM_ARRAY_TASK_ID - 1) / 3 + 1 ))
CELLTYPE=${CELLS[$INDC]}

GENOMES=( NULL hg38 mm10 rheMac10 )
(( INDG=($SLURM_ARRAY_TASK_ID - 1) % 3 + 1 ))
GENOME=${GENOMES[$INDG]}

MEMEDB=( NULL Homo_sapiens.meme Mus_musculus.meme Macaca_mulatta.meme )
DB=${MEMEDB[$INDG]}

##################################################### 
# generate negatives from hg38, mm10, and rheMac10 ##
FGFASTA=$DATADIR/fasta/${GENOME}_${CELLTYPE}_positive.fa.gz
BGFASTA=$DATADIR/fasta/${GENOME}_${CELLTYPE}_biasAway10x.fa

DBDIR=/home/bnphan/resources/motif_databases/CIS-BP
OUTDIR=${DATADIR}/meme-chip; mkdir -p $OUTDIR

POSFILE=$DATADIR/fasta/${GENOME}_${CELLTYPE}_positive.fa
if [[ ! -f $POSFILE ]]; then zcat $FGFASTA > $POSFILE; fi

/home/ikaplow/anaconda2/bin/meme-chip \
-spamo-skip -fimo-skip -o ${OUTDIR}/${GENOME}.${CELLTYPE} \
-db ${DBDIR}/${DB} ${POSFILE}




