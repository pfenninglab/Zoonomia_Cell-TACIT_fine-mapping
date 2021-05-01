#!/bin/bash
#SBATCH --partition=pool1,pfen1,pfen_bigmem
#SBATCH --time 24:00:00
#SBATCH --job-name=Bsnpvars
#SBATCH --dependency=afterok:1536668
#SBATCH --mem=45G
#SBATCH --error=logs/calc_snpvars_%A_%a.txt
#SBATCH --output=logs/calc_snpvars_%A_%a.txt
#SBATCH --array=1-43%4

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
CUTOFF=5e-8

echo "Working on ${PREFIX} GWAS with P < ${CUTOFF} cutoff for loci."

##############################################################################################
# https://github.com/omerwe/polyfun/wiki/1.-Computing-prior-causal-probabilities-with-PolyFun
SUMSTATS=$DATADIR/munged/${PREFIX}.parquet
OUTDIR=${DATADIR}/${PREFIX}/snpvars
mkdir -p $OUTDIR $DATADIR $CACHEDIR

# 3. Compute LD-scores for each SNP bin, do per chromosome
for CHR in {1..22}; do
if [ ! -f ${OUTDIR}/${PREFIX}.${CHR}.l2.ldscore.parquet ]; then
python ${POLYFUNDIR}/polyfun.py --compute-ldscores --chr ${CHR} \
--output-prefix ${OUTDIR}/$PREFIX --ld-ukb --ld-dir $CACHEDIR
fi
done

# 4. Re-estimate per-SNP heritabilities via S-LDSC
if [ ! -f ${OUTDIR}/${PREFIX}.22.snpvar_constrained.gz ]; then
python ${POLYFUNDIR}/polyfun.py --compute-h2-bins \
--allow-missing --output-prefix ${OUTDIR}/$PREFIX \
--sumstats $SUMSTATS --w-ld-chr ${ANNOTDIR}/weights.UKB.
fi
