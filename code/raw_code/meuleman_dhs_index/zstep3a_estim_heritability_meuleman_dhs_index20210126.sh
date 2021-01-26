#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --time 1-0
#SBATCH --job-name=est_herit
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --error=logs/estim_heritability_%A_%a.err.txt
#SBATCH --output=logs/estim_heritability_%A_%a.out.txt
#SBATCH --array=1-65

log2results() {
	IN=$1
	OUT=$(echo $IN| sed 's/.log$/.results/g')
	awk -F '\t' '/Total Observed scale h2*/{flag=1;next}/Lambda GC/{flag=0}flag' $IN | \
	# append line with leading white space to previous
	sed ':a;N;$!ba;s/\n //g'| \
	# parse Partitioned heritability rownames, starts w/ colon
	awk -F":" -v OFS='\t' '{gsub("[[:space:]]", "_", $1); print}' | \
	# transpose row to columns
	awk -v OFS='\t' '
	{ gsub(":", "")
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
	}' | \
	awk -v OFS='\t' '{gsub("_[0-9]+$", "", $1); gsub("L2$", "", $1); print}' | \
	# Add in Observed Scale heritability Z-score
	awk -v OFS='\t'  '{s=(NR==1)?"Observed_scale_z":$2/$3;$0=$0 OFS s}1' | \
	# Add in Conditional Regression Tau Coefficient Z-score
	awk -v OFS='\t'  '{s=(NR==1)?"Coefficient_z":$7/$8;$0=$0 OFS s}1' | \
	awk -v OFS='\t'  '{print $1, $2, $3, $9, $4, $5, $6, $7, $8, $10 }' > ${OUT}
}

# set up the different directors of where files are
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
SNPLIST=${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/listHM3.noMHC.txt
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/meuleman_dhs_index; DATADIR=${SETDIR}/data/raw_data/meuleman_dhs_index
cd $CODEDIR
source activate ldsc; mkdir -p ${DATADIR}/enrichments

# get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' ${GWASDIR}/gwas_list_sumstats.tsv)
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${GWASDIR}/gwas_list_sumstats.tsv)
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' ${GWASDIR}/gwas_list_sumstats.tsv)
CELLTYPES=$(stat -c "%n" $DATADIR/peaks/*.narrowPeak.gz | sed 's!.*/!!' | cut -d'.' -f2 | sort | uniq | sed '/All/d' )
PEAKTYPES=$(stat -c "%n" $DATADIR/peaks/*.narrowPeak.gz | sed 's!.*/!!' | cut -d'.' -f1 | sort | uniq )

for PEAKTYPE in $PEAKTYPES; do
	## do for the full-length peaks
	BACKGROUND=$(stat -c "%n" $DATADIR/annotation/${PEAKTYPE}.All.full.*${POP}.1.annot.gz | sed "s/${POP}.*/${POP}.,/g" | uniq | sed -z "s/\n//g"  )
	for CELL in $CELLTYPES; do
		## run LD score regression to get estimate SNP heritability 
		FOREGROUND=$(stat -c "%n" $DATADIR/annotation/${PEAKTYPE}.${CELL}.full.*${POP}.1.annot.gz | sed "s/${POP}.*/${POP}.,/g" | uniq | sed -z "s/\n//g"  )
		ldsc.py --h2 $GWAS --print-coefficients --print-snps ${SNPLIST} \
		--w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
		--ref-ld-chr ${FOREGROUND}${BACKGROUND}${GWASDIR}/1000G_ALL_Phase3_hg38_files/baselineLD_v2.2/baselineLD_v2.2.${POP}. \
		--frqfile-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/plink/1000G.${POP}.HM3. \
		--out ${DATADIR}/enrichments/${PEAKTYPE}.${CELL}.${GWAS_Label}.${POP}
		## format logs file to make .results file, --print-snps doesn't make a .results file
		log2results ${DATADIR}/enrichments/${PEAKTYPE}.${CELL}.${GWAS_Label}.${POP}.log
	done

	## do for the core peaks
	BACKGROUND=$(stat -c "%n" $DATADIR/annotation/${PEAKTYPE}.All.core.*${POP}.1.annot.gz | sed "s/${POP}.*/${POP}.,/g" | uniq | sed -z "s/\n//g"  )
	for CELL in $CELLTYPES; do
		## run LD score regression to get estimate SNP heritability 
		FOREGROUND=$(stat -c "%n" $DATADIR/annotation/${PEAKTYPE}.${CELL}.core.*${POP}.1.annot.gz | sed "s/${POP}.*/${POP}.,/g" | uniq | sed -z "s/\n//g"  )
		ldsc.py --h2 $GWAS --print-coefficients --print-snps ${SNPLIST} \
		--w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
		--ref-ld-chr ${FOREGROUND}${BACKGROUND}${GWASDIR}/1000G_ALL_Phase3_hg38_files/baselineLD_v2.2/baselineLD_v2.2.${POP}. \
		--frqfile-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/plink/1000G.${POP}.HM3. \
		--out ${DATADIR}/enrichments/${PEAKTYPE}.${CELL}.${GWAS_Label}.${POP}
		## format logs file to make .results file, --print-snps doesn't make a .results file
		log2results ${DATADIR}/enrichments/${PEAKTYPE}.${CELL}.${GWAS_Label}.${POP}.log
	done
done
