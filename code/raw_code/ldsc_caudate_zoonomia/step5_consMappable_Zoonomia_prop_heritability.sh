#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --time=24:00:00
#SBATCH --job-name=cons_herit
#SBATCH --dependency=afterany:1457431
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --error=logs/est_herit_conserved_mappable_%A_%a.txt
#SBATCH --output=logs/est_herit_conserved_mappable_%A_%a.txt
#SBATCH --array=1-64

log2results() {
awk -F '\t' '/Total Observed scale h2*/{flag=1;next}/Lambda GC/{flag=0}flag' ${OUT}.log | \
# append line with leading white space to previous
sed ':a;N;$!ba;s/\n //g'| \
# parse Partitioned heritability rownames, starts w/ colon
awk -F":" -v OFS='\t' '{gsub("[[:space:]]", "_", $1);gsub(":", "", $0); print}' | \
# transpose row to columns
awk -v OFS='\t' '
{ 
for (i=1; i<=NF; i++){
a[NR,i] = $i
}
}
NF>p { p = NF }
END {    
for(j=1; j<=p; j++){
str=a[1,j]
for(i=2; i<=NR; i++){
str=str"\t"a[i,j];
}
print str
}
}' | awk -v OFS='\t' '{gsub("_[0-9]+$", "", $1); gsub("L2$", "", $1); print}' | \
gzip > ${OUT}.results.gz
}

# get the GWAS for this array job
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=${SETDIR}/data/tidy_data/ldsc_gwas; SNPLIST=${GWASDIR}/listHM3.noMHC.txt
CODEDIR=${SETDIR}/code/raw_code/ldsc_caudate_zoonomia; 
DATADIR=${SETDIR}/data/raw_data/ldsc_caudate_zoonomia
cd $CODEDIR; 
source activate ldsc

# get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' ${GWASDIR}/gwas_list_sumstats.tsv)
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${GWASDIR}/gwas_list_sumstats.tsv)
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' ${GWASDIR}/gwas_list_sumstats.tsv)
OUTDIR=${DATADIR}/prop_herit_phyloP; mkdir -p $OUTDIR

###################################################################
# get the cell types that are signif in caudate for this gwas
CELLS=$(awk -v LOOK=${GWAS_Label} '$1 ~ LOOK {print $2; exit}' ${SETDIR}/data/raw_data/caudate_conservation_ldsc/rdas/caudate_celltype_gwas_enriched.tsv | \
sed 's/"//g;s/\r//g' | awk '{$1=$1;print}' | tr ',' ' ')
CELLS="MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Oligo OPC Microglia"
echo "This GWAS is enriched in ${CELLS}."

###################################################################
# run LD score regression over the Corces caudate binary annotations
for CELL in $CELLS; do
CTS_FN=${SETDIR}/data/tidy_data/Zoonomia_data/caudate_zoonomia_mapped_phyloP.${CELL}.${POP}.hg38.celltypes.ldcts
for IND in $(wc -l $CTS_FN | cut -d ' ' -f1 | xargs seq ); do
## use DHS/roadmap/cCRE/Corces2020 merged peaks as bg instead of All
LAB=$(awk -v IND=$IND 'FNR == IND {print $1}' $CTS_FN)
SPECIES=$(echo $LAB | cut -d '.' -f 4)
if [[ $SPECIES == 'Homo_sapiens' || $SPECIES == 'Macaca_mulatta' || $SPECIES == 'Mus_musculus' || $SPECIES == 'Rousettus_aegyptiacus' ]]; then
CELLTYPES=$(awk -v IND=$IND 'FNR == IND {print $2}' $CTS_FN)
CELL2=$(echo $CELLTYPES | cut -d ',' -f1)
OUT=${OUTDIR}/caudate_zoonomia_mapped_phyloP.${GWAS_Label}.${POP}.${CELL}.${LAB}
if [[ ! -f "${OUT}.results.gz" || "${OUT}.results.gz" -ot ${CELL2}.1.l2.M ]]; then
ldsc.py --h2 $GWAS --w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
--ref-ld-chr ${CELLTYPES},${GWASDIR}/1000G_ALL_Phase3_hg38_files/baselineLD_v2.2/baselineLD_v2.2.${POP}. \
--print-coefficients --out ${OUT} 
log2results; rm ${OUT}.log
fi
fi
done

## extract only the top cell type enrichment
if [[ ! -f ${OUTDIR}/caudate_zoonomia_mapped_phyloP.${GWAS_Label}.${POP}.${CELL}.agg.gz ]]; then
gunzip ${OUTDIR}/caudate_zoonomia_mapped_phyloP.${GWAS_Label}.${POP}.${CELL}*.results.gz
awk ' FNR == 2 || NR==1 {print}' \
	${OUTDIR}/caudate_zoonomia_mapped_phyloP.${GWAS_Label}.${POP}.${CELL}*.results |
	gzip >  ${OUTDIR}/caudate_zoonomia_mapped_phyloP.${GWAS_Label}.${POP}.${CELL}.agg.gz
gzip ${OUTDIR}/caudate_zoonomia_mapped_phyloP.${GWAS_Label}.${POP}.${CELL}*.results
fi
done
