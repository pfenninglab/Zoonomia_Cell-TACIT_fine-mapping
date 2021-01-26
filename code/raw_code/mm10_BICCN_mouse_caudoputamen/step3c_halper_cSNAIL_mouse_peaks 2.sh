SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/mm10_BICCN_mouse_caudoputamen
DATADIR=${SETDIR}/data/raw_data/mm10/Mouse_cSNAIL_D1D2
cd $CODEDIR

#####################
# annotate for LDSC
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
HALDIR=${DATADIR}/halper
mkdir -p $HALDIR

for PEAK in ${DATADIR}/peak/*.narrowPeak.gz; do
# zcat $PEAK | awk -v OFS='\t' '{$4="."; print}'| gzip > tmp.txt.gz
# mv tmp.txt.gz $PEAK
sbatch --mem 10G -p pfen1 -w compute-1-40 \
${SETDIR}/code/raw_code/hal_scripts/halper_map_peak_orthologs.sh \
-b $PEAK -o $HALDIR -s Mus_musculus -t Homo_sapiens
done
