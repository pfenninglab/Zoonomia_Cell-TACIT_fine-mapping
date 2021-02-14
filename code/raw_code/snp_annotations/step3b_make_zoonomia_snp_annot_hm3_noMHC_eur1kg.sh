#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
#SBATCH --job-name=zoo_snps
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --error=logs/zoo_snps_%A_%a.err.txt
#SBATCH --output=logs/zoo_snps_%A_%a.out.txt
#SBATCH --array=1-22

SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
SNPLIST=${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/list.txt
ANNOTDIR=${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/zoonomia_snps
DATADIR=${SETDIR}/data/raw_data/snp_annotations
OUTDIR=${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/zoonomia_snps
CHR=${SLURM_ARRAY_TASK_ID}

cd ${SETDIR}/code/raw_code/snp_annotations
mkdir -p ${ANNOTDIR}; rm -f ${ANNOTDIR}/zooSNP.*.hm3_noMHC.${CHR}.tmp.txt

#### Get the SNPs that need to be extracted #####
EUR_SNP=${ANNOTDIR}/tmp1.${CHR}.EUR_SNP.txt
cat ${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/plink_files/1000G.EUR.hg38.${CHR}.bim |\
	cut -f1-4 | awk -v OFS='\t' '{$5=++c; print $1, $2, $3, $4, $5}' | sort -k2,2 > $EUR_SNP

### make lookup table of SNPs and SNP names
BEDFILE=${SETDIR}/data/raw_data/snp_annotations/bed/UKB_imputed_snps.GRCh38.${CHR}.bed.gz
LOOKUP=${ANNOTDIR}/zooSNP.hm3_noMHC.${CHR}.lookup.txt
cat <(zcat $BEDFILE | awk -v OFS='\t' '{print $1 ":" $2 "-" $3, $4}') \
	<(zcat $BEDFILE | awk -v OFS='\t' '{print $4, $4}') | sort -k1,1 > $LOOKUP

### get the list of species in Zoonomia ##
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
TARGETS=$(awk -F'\t' -v OFS='\t' 'FNR >1 {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv | \
	sed '/Homo_sapiens/d' | sed '/Carlito_syrichta/d' | sed '/Manis_tricuspis/d' | sort)

for SPECIES in $TARGETS; do
	## read in the SNPs and get names (sometimes chr:start-end, or sometimes rsID)
	echo "Extracting mapped SNPs on chr${CHR} from ${SPECIES}."
	MAP1=${DATADIR}/halper_zoo/UKB_imputed_snps.GRCh38.${CHR}.Homo_sapiensTo${SPECIES}.halLiftover.tFile.bed.gz
	cat <(zcat $MAP1| cut -f4 | sed 's/:.*//g' ) <(zcat $MAP1| cut -f4 ) | sort | uniq -u  > ${ANNOTDIR}/tmp1.${CHR}.${SPECIES}.txt

	## use look up table to extract rsID's
	join -j 1 -o 1.2 $LOOKUP ${ANNOTDIR}/tmp1.${CHR}.${SPECIES}.txt | \
		sed 's/_.*//g' |  sort | uniq -u | awk -v OFS='\t' '{$2=1; print $1, $2}' \
		> ${ANNOTDIR}/tmp2.${CHR}.${SPECIES}.txt

	# annotate which HM3 SNPs has been halLiftover-ed in SPECIES
	OUTFN=${ANNOTDIR}/zooSNP.${SPECIES}.hm3_noMHC.${CHR}.tmp.txt
	join -1 2 -a 1 -o 1.5 2.2 -e 0 $EUR_SNP \
	${ANNOTDIR}/tmp2.${CHR}.${SPECIES}.txt | awk '!seen[$1]++' | sort -k1,1n | \
	cut -d' ' -f2 | awk -v OFS='\t' -v VAR=$SPECIES 'BEGIN{print VAR}1' > ${OUTFN}
	rm ${ANNOTDIR}/tmp1.${CHR}.${SPECIES}.txt ${ANNOTDIR}/tmp2.${CHR}.${SPECIES}.txt
done

#######################################################################
## merge the annotations from each species together and add SNP info ##
ANNOTFN=${ANNOTDIR}/Zoonomia_240v2.hm3_noMHC.${CHR}.annot.gz
SPECIES_FN=$(ls ${ANNOTDIR}/zooSNP.*.hm3_noMHC.${CHR}.tmp.txt )
FIRST_COL=${ANNOTDIR}/tmp2.${CHR}.EUR_SNP.txt
sort -k5,5n  $EUR_SNP| cut -f1-4 | sed -e '1i\CHR\tSNP\tCM\tBP' > $FIRST_COL
paste $FIRST_COL $SPECIES_FN |gzip > ${ANNOTFN}

### generate LD scores for Zoonomia SNPs 
source activate ldsc
~/src/ldsc/ldsc.py --l2 --ld-wind-kb 1000 --print-snps ${SNPLIST} \
--bfile ${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/plink_files/1000G.EUR.hg38.${CHR} \
--out ${ANNOTDIR}/Zoonomia_240v2.hm3_noMHC.${CHR} --annot ${ANNOTFN}

### last clean up
rm -f ${ANNOTDIR}/zooSNP.*.hm3_noMHC.${CHR}.tmp.txt ${ANNOTDIR}/tmp1.${CHR}.* ${ANNOTDIR}/tmp2.${CHR}.*


