#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --job-name=mk_trackfile
#SBATCH --time 12:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --error=logs/make_trackdb_%A.txt
#SBATCH --output=logs/make_trackdb_%A.txt

### Set up the directories
PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${PROJDIR}/code/final_code/track_files
DATADIR=${PROJDIR}/data/tidy_data/track_files
cd $CODEDIR; mkdir -p $DATADIR/trackdb

source ~/.bashrc


#################################################################
# get the species and cell types that need to be scored by CNNs
CELLS="MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Oligo OPC Microglia"

for CELL in $CELLS; do
echo $CELL
TRACKDB_FN=$DATADIR/trackdb/trackdb_source_${CELL}.txt
OUT="track ${CELL}@@type bigNarrowPeak@@visibility full@@signalFilter 0@@signalFilterLimits 0:10000@@pValueFilter 0@@pValueFilterLimits 0:10000qValueFilter 0@@qValueFilterLimits 0:10000@@shortLabel hs${CELL} bigNPk@@longLabel bigNarrowPeak HomoSapiens_Caudate_${CELL}@@bigDataUrl http://genome.ucsc.edu/goldenPath/help/examples/bigNarrowPeak.bb"
echo $OUT|sed 's/@@/\n/g' > $TRACKDB_FN

TRACKDB_FN=$DATADIR/trackdb/trackdb_ortholog_${CELL}.txt
OUT2="track ${CELL}OrthologPredictionTrack@@type bigNarrowPeak@@visibility full@@signalFilter 0@@signalFilterLimits 0:1@@shortLabel hs${CELL}_orthPred bigNPk@@longLabel bigNarrowPeak HomoSapiens_Caudate_${CELL}_OrthologPredictions@@bigDataUrl http://genome.ucsc.edu/goldenPath/help/examples/bigNarrowPeak.bb"
echo $OUT2|sed 's/@@/\n/g' > $TRACKDB_FN
done