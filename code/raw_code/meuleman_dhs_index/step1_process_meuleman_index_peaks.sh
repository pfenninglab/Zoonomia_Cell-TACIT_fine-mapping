SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/meuleman_dhs_index
DATADIR=${SETDIR}/data/raw_data/meuleman_dhs_index
mv $CODEDIR/annotate* $CODEDIR/logs; mv $DATADIR/annotate* $DATADIR/logs; cd $CODEDIR

#################################################################
# download the Meuleman et al. 2020 hg38 human index DHS regions
mkdir -p $DATADIR/peaks $CODEDIR/logs
# wget https://www.meuleman.org/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz -P $DATADIR/peaks
IDX_PEAKS=${DATADIR}/peaks/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz

# extract list of tissues, change tissue name to lower
zcat $IDX_PEAKS | awk -F'\t' 'OFS= "\t" {$10 = tolower($10); print}' | gzip > tmp.txt.gz
mv tmp.txt.gz $IDX_PEAKS
zcat $IDX_PEAKS | cut -f10 | sed 's/ \/ /\n/g' | sort | uniq | \
	awk -F '\t' 'OFS = "\t" {
		gsub("\\.", "", $1);  
		$2 = $1; 
		gsub(" ", "_", $2); 
	 	gsub("_$", "", $2); print $1, $2}' \
	> ${DATADIR}/DHS_tissue.txt
awk '$0 !~ /component/' ${DATADIR}/DHS_tissue.txt > tmp.txt
mv tmp.txt ${DATADIR}/DHS_tissue.txt

#################################
# split full index DHS by tissue
NUM_TISSUE=$(wc -l ${DATADIR}/DHS_tissue.txt | cut -d ' ' -f1)
for i in `seq $NUM_TISSUE` ; do
	MATCH=$(awk -v IND="$i" 'NR==IND{print $1}' ${DATADIR}/DHS_tissue.txt)
	TISSUE=$(awk -F'\t' -v IND="$i" 'NR==IND{print $2}' ${DATADIR}/DHS_tissue.txt)

	# echo "Extracting full peaks for $TISSUE."
	# ## extract the narrowpeak format for full peak
	# zcat $IDX_PEAKS | grep $MATCH | awk -F '\t' -v OFS='\t' \
	# '{ gsub(" |\\.", "", $10); gsub("/", ",", $10); 
	# print $1, $2, $3, $1":"$2"-"$3":"$10, 0, ".", $5, -1, -1, $7 - $2}' | \
	# gzip > ${DATADIR}/peaks/DHS_Index_and_Vocabulary.${TISSUE}.full.narrowPeak.gz

	echo -e "Extracting core peaks for $TISSUE.\n"
	## extract the narrowpeak format for core peak
	zcat $IDX_PEAKS | grep $MATCH | awk -F '\t' -v OFS='\t' \
	'$8 != $9 && $8 != $7 && $9 != $7 { gsub(" |\\.", "", $10); gsub("/", ",", $10); 
	print $1, $8, $9, $1":"$2"-"$3":"$10, 0, ".", $5, -1, -1, $7 - $8}' | \
	gzip > ${DATADIR}/peaks/DHS_Index_and_Vocabulary.${TISSUE}.core.narrowPeak.gz
done

############################################
## convert all peaks from table to narrowpeak
zcat $IDX_PEAKS | awk -F '\t' -v OFS='\t' 'NR>1 { gsub(" |\\.", "", $10); gsub("/", ",", $10); 
print $1, $2, $3, $1":"$2"-"$3":"$10, 0, ".", $5, -1, -1, $7 - $2}' | \
gzip > ${DATADIR}/peaks/DHS_Index_and_Vocabulary.All.full.narrowPeak.gz

## extract the narrowpeak format for core peak
zcat $IDX_PEAKS | awk -F '\t' -v OFS='\t' \
'NR>1 && $8 != $9 && $8 != $7 && $9 != $7 { gsub(" |\\.", "", $10); gsub("/", ",", $10); 
print $1, $8, $9, $1":"$2"-"$3":"$10, 0, ".", $5, -1, -1, $7 - $8}' | \
gzip > ${DATADIR}/peaks/DHS_Index_and_Vocabulary.All.core.narrowPeak.gz

########################################
# lift halper peaks to rheMac8 and mm10 
HALDIR=${DATADIR}/halper; mkdir -p $HALDIR
for PEAK in ${DATADIR}/peaks/*.narrowPeak.gz; do
	OUTFILE=$(basename $PEAK | sed 's/.narrowPeak.gz/.Homo_sapiensToMus_musculus.HALPER.narrowPeak.gz/g')
	if [ ! -f $HALDIR/${OUTFILE} ]; then
		sbatch --mem 10G -p pfen1 -w compute-1-40 ${SETDIR}/code/raw_code/hal_scripts/halper_map_peak_orthologs.sh \
		-b $PEAK -o $HALDIR -s Homo_sapiens -t Mus_musculus
	fi
done

#####################
# annotate for LDSC
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
ANNOTDIR=${DATADIR}/annotation
mkdir -p $ANNOTDIR

for BED in ${DATADIR}/peaks/*.core.narrowPeak.gz; do
	NAME=$(basename $BED | sed 's/.narrowPeak.gz//g'| sed 's/DHS_Index_and_Vocabulary.//g')
	sbatch --mem 45G --job-name=${NAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
			-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
	sbatch --mem 45G --job-name=${NAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
			-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
done





