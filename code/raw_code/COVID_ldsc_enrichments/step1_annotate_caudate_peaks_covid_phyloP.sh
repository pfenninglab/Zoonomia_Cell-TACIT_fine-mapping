#################
## directories ##
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/COVID_ldsc_enrichments
DATADIR=${SETDIR}/data/raw_data/COVID_ldsc_enrichments
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
PHYLODIR=${GWASDIR}/phyloP/phyloP_cutoff_regions
mkdir -p $DATADIR/annotations $DATADIR/prop_heritability

cd $CODEDIR; mv annot* logs

checkFile(){
HASFILE=TRUE
for CHR in {1..22}; do 
F1=$(echo $FILE |sed "s/#/$CHR/g")
# if no annotation found for chr
if [[ ! -f $F1 ]]; then HASFILE=FALSE
# if annotation crashed and 0 SNPs annotated
elif [[ $(cat $F1) == "0" ]]; then HASFILE=FALSE; fi
done
}

########################################################
# annotate the PhyloP conserved and accelerating regions
for LAB in accl.qval0.05 accl.qval0.01 cons.qval0.05 cons.qval0.01; do
NAME="phyloP.am.${LAB}"
BED=${PHYLODIR}/200m_scoresPhyloP_20210214.${LAB}.bed.gz
FILE=$DATADIR/annotations/${NAME}.#.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
sbatch --mem 24G -p short1,interactive,gpu --time 2:00:00 --job-name=${NAME} \
${CODEDIR}/annotate_for_est_herit_hg38.sh -i ${BED} -n ${NAME} -o $DATADIR/annotations
fi
done





