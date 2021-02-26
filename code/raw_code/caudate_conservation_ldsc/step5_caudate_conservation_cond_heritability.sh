#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen_bigmem
#SBATCH --time=0-8:00:00
#SBATCH --job-name=cond_herit
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --error=logs/cond_heritability_%A_%a.txt
#SBATCH --output=logs/cond_heritability_%A_%a.txt
#SBATCH --array=1-71%10

log2results() {
	awk -F '\t' '/Total Observed scale h2*/{flag=1;next}/Lambda GC/{flag=0}flag' ${OUT}.log | \
	# append line with leading white space to previous
	sed ':a;N;$!ba;s/\n //g'| \
	# parse Partitioned heritability rownames, starts w/ colon
	awk -F":" -v OFS='\t' '{gsub("[[:space:]]", "_", $1);gsub(":", "", $0); print}' | \
	# transpose row to columns
	awk -v OFS='\t' '
	{ 
	    for (i=1; i<=NF; i++)  {
	        a[NR,i] = $i
	    }
	}
	NF>p { p = NF }
	END {    
	    for(j=1; j<=p; j++) {
	        str=a[1,j]
	        for(i=2; i<=NR; i++){
	            str=str"\t"a[i,j];
	        }
	        print str
	    }
	}' | awk -v OFS='\t' '{gsub("_[0-9]+$", "", $1); gsub("L2$", "", $1); print}' | \
	gzip > ${OUT}.results.gz
}

## set up the different directors of where files are
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=${SETDIR}/data/tidy_data/ldsc_gwas; SNPLIST=${GWASDIR}/listHM3.noMHC.txt
CODEDIR=${SETDIR}/code/raw_code/caudate_conservation_ldsc; DATADIR=${SETDIR}/data/raw_data/caudate_conservation_ldsc
OUTDIR=${DATADIR}/conditional_herit
cd $CODEDIR
source activate ldsc; mkdir -p $OUTDIR

## get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' ${GWASDIR}/gwas_list_sumstats.tsv)
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${GWASDIR}/gwas_list_sumstats.tsv)
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' ${GWASDIR}/gwas_list_sumstats.tsv)


#######################################################
## group the peak types e.g. human, human2mouse, etc ##
CTS_FN1=${DATADIR}/caudate_conservation_binary_${POP}_hg38_celltypes.ldcts
TMP_CTS=$( mktemp ${DATADIR}/caudate_conservation_binary_${POP}_hg38_celltypes.XXXXXX )
cat $CTS_FN1 | sed '/^BICCN_CP.MSN_D/d' | sed '/^BICCN_CP.INT/d' | sed '/Consensus/d' | \
	sed '/^200m/d' | sed 's/^BICCN_CP/Mouse_caudate_mappedToHg38/g;s/^Pfenning_Cpu/Mouse_caudate_mappedToHg38/g' |
	sed 's/^Stauffer_caudate/Rhesus_caudate_mappedToHg38/g;s/^Corces2020_caudate/Human_caudate/g' > $TMP_CTS
TYPES=$( cut -f1 $TMP_CTS | cut -d '.' -f1 | sort | uniq )
CELLS=$( cut -f1 $TMP_CTS | cut -d '.' -f2 | sort | uniq )
BGPEAKS=$( cut -f2 $TMP_CTS | cut -d ',' -f2 | sort | uniq )

########################################################
## Conditional SNP heritability on caudate cell types ##
for LAB in $TYPES; do
OUT=${OUTDIR}/caudate_conservation_binary.${LAB}.${GWAS_Label}.${POP}

## gather the foreground cell types for conditional cell type enrichments
CELLTYPES=""
for CELL in $CELLS; do CELLTYPES=${CELLTYPES},$(awk -v VAR="${LAB}.${CELL}" '{if(match($1, VAR)) print $2}' $TMP_CTS | cut -d ',' -f1 ) ; done
CELLTYPES=$(echo $CELLTYPES | sed 's/^,//g')

## perform the conditional cell type specific enrichments
if [[ ! -f "${OUT}.results.gz" ]]; then # --print-snps ${SNPLIST} 
ldsc.py --h2 $GWAS --print-coefficients \
--w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
--ref-ld-chr ${CELLTYPES},${BGPEAKS},${GWASDIR}/1000G_ALL_Phase3_hg38_files/baselineLD_v2.2/baselineLD_v2.2.${POP}. \
--out ${OUT}
log2results
fi
done
rm $TMP_CTS

