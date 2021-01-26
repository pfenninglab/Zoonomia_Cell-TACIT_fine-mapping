#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
##SBATCH --time 1-0
#SBATCH --mem=15G
#SBATCH --job-name=meme
#SBATCH --error=logs/BICCN_mouse_caudoputamen_meme-chip_%A_%a_out.txt
#SBATCH --output=logs/BICCN_mouse_caudoputamen_meme-chip_%A_%a_out.txt
#SBATCH --array=1-10

# cd /projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/code/raw_code/meme-chip_analyses

PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/raw_data/mm10/BICCN_mouse_caudoputamen
GENOME=/home/bnphan/src/atac-seq-pipeline/genome/mm10/mm10_no_alt_analysis_set_ENCODE.fasta
DBFILE=/home/bnphan/resources/motif_databases/MOUSE/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme

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

# make meme-chip output dir and run
OUTDIR=${PROJDIR}/meme/$(echo $NARROWPEAK | sed 's/.narrowPeak.gz//g' )
/home/ikaplow/anaconda2/bin/meme-chip \
	-fimo-skip -spamo-skip -noecho \
	-oc $OUTDIR -db $DBFILE \
	$POSFILE

rm ${TMPBED} ${POSFILE} ${OUTDIR}/*.fasta