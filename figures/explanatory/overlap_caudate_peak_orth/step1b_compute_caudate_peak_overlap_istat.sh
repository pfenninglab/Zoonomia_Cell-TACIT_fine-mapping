SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/figures/explanatory/overlap_caudate_peak_orth
DATADIR=${SETDIR}/data/raw_data/hg38/Corces_2020
MACDIR=${SETDIR}/data/raw_data/rheMac10/Stauffer_caudate
MUSDIR=${SETDIR}/data/raw_data/mm10/BICCN_mouse_caudoputamen

cd ${CODEDIR}; mkdir -p ${CODEDIR}/logs ${CODEDIR}/tables

for REF in ${DATADIR}/peak/Corces2020_caudate.*.narrowPeak.gz; do
	for QUERY in ${DATADIR}/peak/*.narrowPeak.gz ${MUSDIR}/halper/*ToHomo_sapiens.HALPER.narrowPeak.gz ${MACDIR}/halper/*ToHomo_sapiens.HALPER.narrowPeak.gz; do
		echo "Reference: ${REF}		Query: ${REF}"
		sbatch --mem 3G -p pool1 --time 1-0 ${CODEDIR}/compute_overlap_istat.sh -q ${QUERY} -r ${REF} -o ${CODEDIR}/tables
	done
done
