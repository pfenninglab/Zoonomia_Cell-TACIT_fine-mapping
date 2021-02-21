#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=interactive
#SBATCH --time 0-5
#SBATCH --job-name=zoo_snps
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --error=logs/zoo_snps_%A_%a.err.txt
#SBATCH --output=logs/zoo_snps_%A_%a.out.txt
#SBATCH --array=1-21

SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=${SETDIR}/data/tidy_data/ldsc_gwas
SNPLIST=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/listHM3.noMHC.txt
ANNOTDIR=${GWASDIR}/1000G_ALL_Phase3_hg38_files/zoonomia_snps
DATADIR=${SETDIR}/data/raw_data/snp_annotations
OUTDIR=${GWASDIR}/1000G_ALL_Phase3_hg38_files/zoonomia_snps
CHR=${SLURM_ARRAY_TASK_ID}

cd ${SETDIR}/code/raw_code/snp_annotations
mkdir -p ${ANNOTDIR} 

#### Get the SNPs that need to be extracted #####
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
AFR_SNP=${ANNOTDIR}/tmp1.${CHR}.AFR_SNP.txt; EUR_SNP=${ANNOTDIR}/tmp1.${CHR}.EUR_SNP.txt
zcat ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.AFR.weights.hm3_noMHC.${CHR}.l2.ldscore.gz | cut -f1-3 | awk -v OFS='\t' 'NR>1{print}' | sort -k2,2 > $AFR_SNP
zcat ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.EUR.weights.hm3_noMHC.${CHR}.l2.ldscore.gz | cut -f1-3 | awk -v OFS='\t' 'NR>1{print}' | sort -k2,2 > $EUR_SNP

### make lookup table of SNPs and SNP names
BEDFILE=${SETDIR}/data/raw_data/snp_annotations/bed/UKB_imputed_snps.GRCh38.${CHR}.bed.gz
LOOKUP=${ANNOTDIR}/zooSNP.SPECIES.hm3_noMHC.${POP}.${CHR}.lookup.txt
cat <(zcat $BEDFILE | awk -v OFS='\t' '{print $1 ":" $2 "-" $3, $4}') \
	<(zcat $BEDFILE | awk -v OFS='\t' '{print $4, $4}') | sort -k1,1 > $LOOKUP

### get the list of species in Zoonomia ##
TARGETS=$(awk -F'\t' -v OFS='\t' 'FNR >1 {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv | \
	sed '/Homo_sapiens/d' | sed '/Carlito_syrichta/d' | sed '/Manis_tricuspis/d' | sort)
## clear up any old runs
rm -f ${ANNOTDIR}/zooSNP.*.hm3_noMHC.${POP}.${CHR}.tmp.txt

for SPECIES in $TARGETS; do
## read in the SNPs and get names (sometimes chr:start-end, or sometimes rsID)
echo "Extracting mapped SNPs on chr${CHR} from ${SPECIES}."
MAP1=${DATADIR}/halper_zoo/UKB_imputed_snps.GRCh38.${CHR}.Homo_sapiensTo${SPECIES}.halLiftover.tFile.bed.gz
cat <(zcat $MAP1| cut -f4 | sed 's/:.*//g' ) <(zcat $MAP1| cut -f4 ) | uniq -u | sort \
	> ${ANNOTDIR}/tmp1.${CHR}.${SPECIES}.txt

## use look up table to extract rsID's
join -j 1 -o 1.2 $LOOKUP ${ANNOTDIR}/tmp1.${CHR}.${SPECIES}.txt | \
	sed 's/_.*//g' | uniq -u | sort | awk -v OFS='\t' '{$2=1; print $1, $2}' \
	> ${ANNOTDIR}/tmp2.${CHR}.${SPECIES}.txt

for POP in AFR EUR; do
	# annotate which HM3 SNPs has been halLiftover-ed in SPECIES
	OUTFN=${ANNOTDIR}/zooSNP.${SPECIES}.hm3_noMHC.${POP}.${CHR}.tmp.txt
	join -1 2 -a 1 -o 1.1 1.2 1.3 2.2 -e 0 ${ANNOTDIR}/tmp1.${CHR}.${POP}_SNP.txt \
	${ANNOTDIR}/tmp2.${CHR}.${SPECIES}.txt | awk '!seen[$2]++' | sort -k1,1 -k3,3n | \
	cut -d' ' -f4 | awk -v OFS='\t' -v VAR=$SPECIES 'BEGIN{print VAR}1' > ${OUTFN}
done
rm ${ANNOTDIR}/tmp1.${CHR}.${SPECIES}.txt ${ANNOTDIR}/tmp2.${CHR}.${SPECIES}.txt
done

#######################################################################
## merge the annotations from each species together and add SNP info ##
source activate ldsc
for POP in AFR EUR; do
ANNOTFN=${ANNOTDIR}/Zoonomia_240v2.hm3_noMHC.${POP}.${CHR}.annot.gz
SPECIES_FN=$(ls ${ANNOTDIR}/zooSNP.*.hm3_noMHC.${POP}.${CHR}.tmp.txt )
FIRST_COL=${ANNOTDIR}/tmp2.${CHR}.EUR_SNP.txt
awk '!seen[$2]++ {$4=0; print $1, $2, $3, $4}' ${ANNOTDIR}/tmp1.${CHR}.${POP}_SNP.txt | \
sort -k1,1 -k3,3n | sed -e '1i\CHR\tSNP\tBP\tCM' > $FIRST_COL
paste $FIRST_COL $SPECIES_FN |gzip > ${ANNOTFN}

### generate LD scores for Zoonomia SNPs
~/src/ldsc/ldsc.py --l2 --print-snps ${SNPLIST} --ld-wind-kb 1000 \
--bfile /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_ALL_Phase3_hg38_files/plink/1000G.${POP}.HM3.${CHR} \
--out ${ANNOTDIR}/Zoonomia_240v2.hm3_noMHC.v2.${POP}.${CHR} --annot ${ANNOTFN}
done

## last clean up
rm -f ${ANNOTDIR}/zooSNP.*.hm3_noMHC.${POP}.${CHR}.tmp.txt ${ANNOTDIR}/tmp1.${CHR}.* ${ANNOTDIR}/tmp2.${CHR}.*


