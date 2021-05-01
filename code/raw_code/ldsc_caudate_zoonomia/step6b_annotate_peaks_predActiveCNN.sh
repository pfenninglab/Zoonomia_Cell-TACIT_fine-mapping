#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
#SBATCH --job-name=annotScored
#SBATCH --time 12:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --error=logs/score_mappable_peaks_%A_%a.txt
#SBATCH --output=logs/score_mappable_peaks_%A_%a.txt
#SBATCH --array=1-240

checkFile(){
HASFILE=TRUE
for CHR in {1..22}; do 
F1=$(echo $FILE |sed "s/#/$CHR/g")
if [[ ! -f $F1 ]]; then 
# if no annotation found for chr
HASFILE=FALSE
elif [[ $(cat $F1) == "0" ]]; then
# if annotation crashed and 0 SNPs annotated
HASFILE=FALSE
elif [[ $F1 -ot $BED ]]; then
# if newer foreground file
HASFILE=FALSE
fi
done
}

### Set up the directories
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
SNPLIST=${GWASDIR}/listHM3.noMHC.txt
CODEDIR=${SETDIR}/code/raw_code/ldsc_caudate_zoonomia; 
DATADIR=${SETDIR}/data/raw_data/ldsc_caudate_zoonomia
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
PEAKDIR=${SETDIR}/data/raw_data/hg38/Corces_2020/halper_zoo
ANNOTDIR=$DATADIR/annotations
cd $CODEDIR; mkdir -p $ANNOTDIR
source ~/.bashrc

for SLURM_ARRAY_TASK_ID in {1..240}; do

#################################################################
# get the species and cell types that need to be scored by CNNs
CELLS="MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Oligo OPC Microglia"
SPECIES=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv )

##############################################
# Predict SPECIES ortholog CELLTYPE activity
for CELL in $CELLS; do
	echo "Annotating peaks mapped to ${SPECIES} for ${CELL}."
	# 5) get the hg38 coorindates of predicted active peaks in this species
	BED=${DATADIR}/peaks/Corces2020_caudate.${CELL}.Homo_sapiensTo${SPECIES}.predActive.bed.gz
	NAME=Corces2020_caudate.${CELL}.PredActiveIn${SPECIES}

	if [[ -f $BED ]]; then
		## for AFR annotations
		FILE=${ANNOTDIR}/${NAME}.AFR.#.l2.M; checkFile
		if [[ $HASFILE == "FALSE" ]]; then 
		sbatch --mem 4G -p short1,interactive,pool1,pool3-bigmem,pfen1 --array=1-22 \
		--time 2:00:00 --output=/dev/null --error=/dev/null --job-name=AFR.${NAME} \
		${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
		-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
		# sleep 1m
		else echo "Annotations exist for ${NAME} in AFR populations."
		fi

		# for EUR annotations
		FILE=${ANNOTDIR}/${NAME}.EUR.#.l2.M; checkFile
		if [[ $HASFILE == "FALSE" ]]; then 
		sbatch --mem 4G -p short1,interactive,pool1,pool3-bigmem,pfen1 --array=1-22 \
		--time 2:00:00 --output=/dev/null --error=/dev/null --job-name=EUR.${NAME} \
		${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
		-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
		# sleep 1m
		else echo "Annotations exist for ${NAME} in EUR populations."; 
		fi
	else echo "The file ${BED} does not exist."; 
	fi
done

done


