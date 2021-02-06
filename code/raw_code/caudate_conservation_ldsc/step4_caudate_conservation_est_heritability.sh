#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --time=1-0
#SBATCH --job-name=est_herit
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --error=logs/estim_heritability_%A_%a.txt
#SBATCH --output=logs/estim_heritability_%A_%a.txt
#SBATCH --array=48-71

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
cd $CODEDIR
source activate ldsc; mkdir -p ${DATADIR}/est_herit_cond

## get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' ${GWASDIR}/gwas_list_sumstats.tsv)
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${GWASDIR}/gwas_list_sumstats.tsv)
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' ${GWASDIR}/gwas_list_sumstats.tsv)
PEAKTYPES=$(stat -c "%n" $DATADIR/annotation/*.1.annot.gz |grep -E 'Corces2020'| sed 's!.*/!!' | cut -d'.' -f1 | sort | uniq )

########################################################################
## Estimate SNP heritability on binary Meuleman DHS index annotations ##
for PEAKTYPE in $PEAKTYPES; do
	## use DHS/roadmap/cCRE/Corces2020 merged peaks as bg instead of All
	BINARY_DHS_BG=${SETDIR}/data/raw_data/caudate_conservation_ldsc/annotation/Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer.BG.${POP}.,
	CELLTYPES=$(stat -c "%n" $DATADIR/annotation/${PEAKTYPE}.*.${POP}.1.annot.gz | sed "s/${POP}.*/${POP}.,/g" | sed '/Consensus/d' | sed '/All/d' | sort| uniq | sed -z 's/\n//g')
	OUT=${DATADIR}/est_herit_cond/${PEAKTYPE}.binary.${GWAS_Label}.${POP}
	# if [[ ! -f "${OUT}.results.gz" ]]; then
		ldsc.py --h2 $GWAS --print-coefficients --print-snps ${SNPLIST} --n-blocks 400 \
		--w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
		--ref-ld-chr ${CELLTYPES}${BINARY_DHS_BG}${GWASDIR}/1000G_ALL_Phase3_hg38_files/baselineLD_v2.2/baselineLD_v2.2.${POP}. \
		--frqfile-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/plink/1000G.${POP}.HM3. \
		--out ${OUT} 
	# fi
	log2results
	## format logs file to make .results file, --print-snps doesn't make a .results file

	# #############################################################################################
	# ## Estimate SNP heritability on continuous phyloP scores intersected w/ Meuleman DHS peaks ##
	# CONTIN_DHS_BG=${SETDIR}/data/raw_data/caudate_conservation_ldsc/annot_phyloP/Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer.BG.${POP}.,
	# CELLTYPES=$(stat -c "%n" $DATADIR/annot_phyloP/${PEAKTYPE}.*.${POP}.1.annot.gz | sed "s/${POP}.*/${POP}.,/g" | sed '/Consensus/d'| sed '/All/d' | sort| uniq | sed -z 's/\n//g')
	# OUT=${DATADIR}/est_herit_cond/${PEAKTYPE}.phyloP.${GWAS_Label}.${POP}
	# if [[ ! -f "${OUT}.results.gz" ]]; then
	# 	ldsc.py --h2 $GWAS --print-coefficients --print-snps ${SNPLIST} \
	# 	--w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
	# 	--ref-ld-chr ${CELLTYPES}${CONTIN_DHS_BG}${GWASDIR}/1000G_ALL_Phase3_hg38_files/baselineLD_v2.2/baselineLD_v2.2.${POP}. \
	# 	--frqfile-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/plink/1000G.${POP}.HM3. \
	# 	--out ${OUT}  
	# fi
	# log2results
	## format logs file to make .results file, --print-snps doesn't make a .results file
done


