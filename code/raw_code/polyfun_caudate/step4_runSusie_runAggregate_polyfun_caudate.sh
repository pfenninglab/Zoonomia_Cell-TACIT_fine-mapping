
SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate'
CACHEDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/polyfun/LD_cache
ANNOTDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/polyfun/baselineLF2.2.UKB
ZOONOMIADIR=/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping/data/raw_data/zoonomia_annotations/annotation
DATADIR=${SETWD}/data/raw_data/polyfun_caudate; CODEDIR=${SETWD}/code/raw_code/polyfun_caudate
POLYFUNDIR='/home/bnphan/src/polyfun'
cd $CODEDIR;


######################################
## split finemapping job by chromosome
for ID in {1..43}; do
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $3}' ${SETWD}/data/tidy_data/ldsc_gwas/gwas_list_sumstats_polyfun.tsv)
N=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $10}' ${SETWD}/data/tidy_data/ldsc_gwas/gwas_list_sumstats_polyfun.tsv)
CUTOFF=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $11}' ${SETWD}/data/tidy_data/ldsc_gwas/gwas_list_sumstats_polyfun.tsv); 
OUTDIR=${DATADIR}/${PREFIX}/susie 
# 23G good for blocks w/ < 15k SNPs, 47G good for blocks <30k snps
if [[ ! -f ${DATADIR}/${PREFIX}/${PREFIX}.caudate_conservation.aggregate.txt.gz && -f ${OUTDIR}/polyfun_all_jobs_22.txt ]]; then
echo "Submitting job for ${PREFIX} GWAS."
sbatch -p pool1,interactive,pool3-bigmem,gpu,pfen1,pfen_bigmem --mem 23G --time 8:00:00 --array 1-22 --export=JobsFileName="${OUTDIR}/polyfun_all_jobs_@.txt" ${CODEDIR}/slurm_finemap_byLine.sh
fi; done


##########################################
## combine all the jobs together per trait
## catch the jobs that need a lot of RAM
for ID in {1..43}; do
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $3}' ${SETWD}/data/tidy_data/ldsc_gwas/gwas_list_sumstats_polyfun.tsv)
N=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $10}' ${SETWD}/data/tidy_data/ldsc_gwas/gwas_list_sumstats_polyfun.tsv)
OUTDIR=${DATADIR}/${PREFIX}/susie; JobsFileName=${OUTDIR}/polyfun_all_jobs.txt
## merge all jobs together
if [[ ! -f ${DATADIR}/${PREFIX}/${PREFIX}.caudate_conservation.aggregate.txt.gz && -f ${OUTDIR}/polyfun_all_jobs_22.txt ]]; then cat ${OUTDIR}/polyfun_all_jobs_*.txt > $JobsFileName; fi
if [[ ! -f ${DATADIR}/${PREFIX}/${PREFIX}.caudate_conservation.aggregate.txt.gz && -f ${JobsFileName} ]]; then
echo "Submitting job for ${PREFIX} GWAS."; 
sbatch -p pfen_bigmem,pool3-bigmem,pfen3 --mem 120G --time 24:00:00 \
--export=JobsFileName=${JobsFileName} ${CODEDIR}/slurm_finemap_byLine.sh
fi; done


##################################################
## submit job to aggregate fine-mapping when done
## catch the jobs that need a lot of RAM
for ID in {1..43}; do
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $3}' ${SETWD}/data/tidy_data/ldsc_gwas/gwas_list_sumstats_polyfun.tsv)
N=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $10}' ${SETWD}/data/tidy_data/ldsc_gwas/gwas_list_sumstats_polyfun.tsv)
CUTOFF=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $11}' ${SETWD}/data/tidy_data/ldsc_gwas/gwas_list_sumstats_polyfun.tsv); 
SUMSTATS=$DATADIR/munged/${PREFIX}.parquet; OUTDIR=${DATADIR}/${PREFIX}/susie
if [[ ! -f ${DATADIR}/${PREFIX}/${PREFIX}.caudate_conservation.causal_set.txt.gz && -f ${OUTDIR}/polyfun_all_jobs_22.txt ]]; then
echo "Aggregating results from ${PREFIX} GWAS with P < ${CUTOFF} cutoff for loci."
sbatch --export=OUTDIR=${OUTDIR},PREFIX=${PREFIX}.caudate_conservation,SUMSTATS=${SUMSTATS},CUTOFF=${CUTOFF} \
--partition pfen3,pfen1 --time 3:00:00 --mem 30G ${CODEDIR}/slurm_polyfun_aggregate.sh
fi; done

