# directories
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/ldsc_atac_and_conservation
DATADIR=${SETDIR}/data/raw_data
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments

cd $CODEDIR

OUTDIR=${DATADIR}/ldsc_atac_and_conservation/peak_meta
mkdir -p $OUTDIR

############################################################################################
# rank the conservation scores by the same number of peaks as peaks mappable and pred active
CONTYPES="91mam_GERP 200mam_PhyloP 43prim_PhastCons"
PEAKTYPES="Enhancer Promoter Other"
CELLTYPES="MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Microglia OPC Oligo"

for CELLTYPE in $CELLTYPES; do
for PEAKTYPE in $PEAKTYPES; do
for BEDFILE in ${SETDIR}/data/raw_data/ldsc_atac_and_conservation/peak/Corces2020.${CELLTYPE}.*.${PEAKTYPE}.predActive.bed.gz; do
NAME=$(basename $BEDFILE .predActive.bed.gz)
NUMPEAK=$(zcat ${BEDFILE} | wc -l )
echo "Working on ${NAME}."
for CONTYPE in $CONTYPES; do
OUTFILE=${OUTDIR}/${NAME}.${CONTYPE}.bed.gz
if [[ ! -f "$OUTFILE" ]]; then 
CONFILE=${DATADIR}/ldsc_atac_and_conservation/peak/Corces2020.${CELLTYPE}.${CONTYPE}_summary.bed
echo "Creating $(basename $OUTFILE)."
## get the top N peaks sort by the 
sort -k5,5nr $CONFILE | head -${NUMPEAK} | gzip > $OUTFILE
fi
done; done; done; done

