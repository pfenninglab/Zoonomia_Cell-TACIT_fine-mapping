#!/bin/bash
#SBATCH --partition=short1,interactive,gpu
#SBATCH --time 2:00:00
#SBATCH --job-name=split_jobs
#SBATCH --mem=10G
#SBATCH --dependency=afterok:1363930
#SBATCH --error=logs/split_jobs_%A_%a.txt
#SBATCH --output=logs/split_jobs_%A_%a.txt
#SBATCH --array=1-43

SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate'
CACHEDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/polyfun/LD_cache
ANNOTDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/polyfun/baselineLF2.2.UKB
ZOONOMIADIR=/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping/data/raw_data/zoonomia_annotations/annotation
DATADIR=${SETWD}/data/raw_data/polyfun_caudate; CODEDIR=${SETWD}/code/raw_code/polyfun_caudate
POLYFUNDIR='/home/bnphan/src/polyfun'

cd $CODEDIR; 
source ~/.bashrc; conda activate polyfun

PREFIX=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == IND + 1 {print $3}' ${SETWD}/data/tidy_data/ldsc_gwas/gwas_list_sumstats_polyfun.tsv)
N=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == IND + 1 {print $10}' ${SETWD}/data/tidy_data/ldsc_gwas/gwas_list_sumstats_polyfun.tsv)
CUTOFF=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == IND + 1 {print $11}' ${SETWD}/data/tidy_data/ldsc_gwas/gwas_list_sumstats_polyfun.tsv); 

echo "Working on ${PREFIX} GWAS with P < ${CUTOFF} cutoff for loci."

################################################################
# https://github.com/omerwe/polyfun/wiki/3.-Functionally-informed-fine-mapping-with-finemapper
# generate scripts to run poly fun mapping at significant loci
OUTDIR=${DATADIR}/${PREFIX}/susie
mkdir -p $OUTDIR $DATADIR $CACHEDIR

for CHR in {1..22}; do
## if jobs not already created
SUMSTATS=${DATADIR}/${PREFIX}/snpvars/${PREFIX}.${CHR}.snpvar_constrained.gz
if [[ ! -f ${OUTDIR}/polyfun_all_jobs_${CHR}.txt && -f ${DATADIR}/${PREFIX}/snpvars/${PREFIX}.${CHR}.snpvar_constrained.gz ]]; then
# if [[ -f ${DATADIR}/${PREFIX}/snpvars/${PREFIX}.${CHR}.snpvar_constrained.gz ]]; then
echo "Working on chr${CHR}."
python ${POLYFUNDIR}/create_finemapper_jobs.py \
--allow-missing --method susie \
--cache-dir ${CACHEDIR} --sumstats ${SUMSTATS} \
--pvalue-cutoff ${CUTOFF} --max-num-causal 10 \
--out-prefix ${OUTDIR}/polyfun_all --n ${N} \
--jobs-file ${OUTDIR}/polyfun_all_jobs_${CHR}.txt \
--memory 5 --chr ${CHR}
fi
done

