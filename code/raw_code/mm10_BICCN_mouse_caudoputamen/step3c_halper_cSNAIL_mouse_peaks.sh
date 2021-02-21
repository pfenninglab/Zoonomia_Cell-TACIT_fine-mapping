SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/mm10_BICCN_mouse_caudoputamen
DATADIR=${SETDIR}/data/raw_data/mm10
cd $CODEDIR

#####################
# annotate for LDSC
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments

for PEAK in ${DATADIR}/Mouse_cSNAIL_D1D2/peak/*narrowPeak.gz; do
	# zcat $PEAK | awk -v OFS='\t' '{$4="."; print}'| gzip > tmp.txt.gz
	# mv tmp.txt.gz $PEAK
	# sbatch --mem 15G -p pool3-bigmem --time 4-4 -w compute-1-35 ${SETDIR}/code/raw_code/hal_scripts/halper_map_peak_orthologs.sh -b $PEAK -o $HALDIR -s Mus_musculus -t Homo_sapiens
	HALDIR=$(dirname $PEAK | sed 's/peak/halper/g'); mkdir -p $HALDIR
	sbatch --mem 8G -p pfen1 --time 4-4 -w compute-1-39 ${SETDIR}/code/raw_code/hal_scripts/halper_map_peak_orthologs.sh \
		-b $PEAK -o $HALDIR -s Mus_musculus -t Macaca_mulatta
	sbatch --mem 8G -p pfen1 --time 4-4 -w compute-1-39 ${SETDIR}/code/raw_code/hal_scripts/halper_map_peak_orthologs.sh \
		-b $PEAK -o $HALDIR -s Mus_musculus -t Homo_sapiens
done
