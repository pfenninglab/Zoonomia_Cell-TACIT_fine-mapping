# directories
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/ldsc_atac_and_conservation
DATADIR=${SETDIR}/data/raw_data
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments

cd $CODEDIR

OUTDIR=${DATADIR}/ldsc_atac_and_conservation/peak
mkdir -p $OUTDIR

###########################################
# get the summary of 91 mammals GERP score
BWFILE=${GWASDIR}/gerp/gerp_conservation_scores.homo_sapiens.GRCh38.bw
CONTYPE=91mam_GERP

for BEDFILE in ${OUTDIR}/Corces2020*main*.bed.gz; do
NAME=$(basename $BEDFILE .bed.gz)
OUTFILE=${OUTDIR}/${NAME}.${CONTYPE}_summary.bed
if [[ ! -f "$OUTFILE" ]]; then 
TMPFILE=${OUTDIR}/$(basename $BEDFILE .gz).txt
## the GERP scores had different chromosome naming system
zcat $BEDFILE| cut -f1-3| sed 's/^chr//g' > $TMPFILE
OUTFILE2=${OUTDIR}/${NAME}.${CONTYPE}_tmp.txt
echo "Working on ${NAME}."
# https://github.com/CRG-Barcelona/bwtool/wiki/summary
bwtool summary $TMPFILE $BWFILE $OUTFILE2 -fill=0
## combine the average conservation scores per chromosome together
## change back the chromosome to UCSC style
cat ${OUTFILE2}| sort -k1,1 -k2,2n | cut -f1-3,8 | \
awk -F'\t' 'BEGIN { OFS = FS} 
{ $1="chr"$1;$5=$2+1;$6="hg38:"$1":"$5"-"$3":250"} 
{ print $1, $2, $3, $6, $4, "."}' \
> ${OUTFILE}
rm $OUTFILE2 $TMPFILE
fi
done



#####################################
# get the summary of 43 primate phastCons
BWFILE=${GWASDIR}/phyloP/Primates_PhastCons_scores_hg38/43prim_PhastCons.bigWig
CONTYPE=43prim_PhastCons

for BEDFILE in ${OUTDIR}/Corces2020*main*.bed.gz; do
NAME=$(basename $BEDFILE .bed.gz)
OUTFILE=${OUTDIR}/${NAME}.${CONTYPE}_summary.bed
if [[ ! -f "$OUTFILE" ]]; then 
TMPFILE=${OUTDIR}/$(basename $BEDFILE .gz).txt
## the GERP scores had different chromosome naming system
zcat $BEDFILE| cut -f1-3| sed 's/^chr//g' > $TMPFILE
OUTFILE2=${OUTDIR}/${NAME}.${CONTYPE}_tmp.txt
echo "Working on ${NAME}."
# https://github.com/CRG-Barcelona/bwtool/wiki/summary
bwtool summary $TMPFILE $BWFILE $OUTFILE2 -fill=0
## combine the average conservation scores per chromosome together
## change back the chromosome to UCSC style
cat ${OUTFILE2}| sort -k1,1 -k2,2n | cut -f1-3,8 | \
awk -F'\t' 'BEGIN { OFS = FS} 
{ $1="chr"$1;$5=$2+1;$6="hg38:"$1":"$5"-"$3":250"} 
{ print $1, $2, $3, $6, $4, "."}' \
> ${OUTFILE}
rm $OUTFILE2 $TMPFILE
fi
done





#####################################
# get the summary of 200mam phyloP
BWFILE=${GWASDIR}/phyloP/human-centered-200m-Feb2021/200m_scoresPhyloP_20210214.bigWig
CONTYPE=200mam_PhyloP

for BEDFILE in ${OUTDIR}/Corces2020*main*.bed.gz; do
NAME=$(basename $BEDFILE .bed.gz)
OUTFILE=${OUTDIR}/${NAME}.${CONTYPE}_summary.bed
if [[ ! -f "$OUTFILE" ]]; then 
TMPFILE=${OUTDIR}/$(basename $BEDFILE .gz).txt
## the GERP scores had different chromosome naming system
zcat $BEDFILE| cut -f1-3| sed 's/^chr//g' > $TMPFILE
OUTFILE2=${OUTDIR}/${NAME}.${CONTYPE}_tmp.txt
echo "Working on ${NAME}."
# https://github.com/CRG-Barcelona/bwtool/wiki/summary
bwtool summary $TMPFILE $BWFILE $OUTFILE2 -fill=0
## combine the average conservation scores per chromosome together
## change back the chromosome to UCSC style
cat ${OUTFILE2}| sort -k1,1 -k2,2n | cut -f1-3,8 | \
awk -F'\t' 'BEGIN { OFS = FS} 
{ $1="chr"$1;$5=$2+1;$6="hg38:"$1":"$5"-"$3":250"} 
{ print $1, $2, $3, $6, $4, "."}' \
> ${OUTFILE}
rm $OUTFILE2 $TMPFILE
fi
done


rm ${OUTDIR}/*tmp*
