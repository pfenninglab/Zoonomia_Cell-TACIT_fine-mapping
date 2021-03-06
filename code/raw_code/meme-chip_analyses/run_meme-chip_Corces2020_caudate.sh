#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
##SBATCH --time 1-0
#SBATCH --mem=15G
#SBATCH --job-name=meme
#SBATCH --error=logs/corces_meme-chip_%A_%a_out.txt
#SBATCH --output=logs/corces_meme-chip_%A_%a_out.txt
#SBATCH --array=1-12

# cd /projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/code/raw_code/meme-chip_analyses

PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/raw_data/hg38/Corces_2020
GENOME=/home/bnphan/resources/genomes/GRCh38.p13/GRCh38.primary_assembly.genome.fa
DBFILE=/home/bnphan/resources/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme

TMPDIR=/scratch/bnphan/meme-chip
mkdir -p $TMPDIR

# get the narrowpeak file
((IND= $SLURM_ARRAY_TASK_ID + 1))
NARROWPEAK=$(ls -l ${PROJDIR}/peak | awk -v IND="${IND}" 'NR==IND {print $9}')
TMPBED=${TMPDIR}/$(echo $NARROWPEAK | sed 's/.narrowPeak.gz/.narrowPeak/g' )
zcat ${PROJDIR}/peak/$NARROWPEAK > ${TMPBED}

# extract fasta sequences
POSFILE=${TMPDIR}/$(echo $NARROWPEAK | sed 's/.narrowPeak.gz/.fasta/g')
bedtools getfasta -name -fi ${GENOME} -bed ${TMPBED} > ${POSFILE} 
rm 

# make meme-chip output dir and run
OUTDIR=${PROJDIR}/meme/$(echo $NARROWPEAK | sed 's/.narrowPeak.gz//g' )
/home/ikaplow/anaconda2/bin/meme-chip \
	-fimo-skip -spamo-skip -noecho \
	-oc $OUTDIR -db $DBFILE \
	$POSFILE

rm ${TMPBED} ${POSFILE} ${OUTDIR}/*.fasta