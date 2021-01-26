function cleanup()
{
    rm ${TMPDIR}/${NAME}*tmp*
}


# default cactus file
# CACTUSFILE=/projects/pfenninggroup/machineLearningForComputationalBiology/alignCactus/10plusway-master.hal
CACTUSFILE=/data/pfenninggroup/machineLearningForComputationalBiology/alignCactus/mam241/241-mammalian-2020v2.hal
OVERWRITE='FALSE'

cd /projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/code/raw_code/meuleman_dhs_index
mv annotate_* logs

# read in command line arguments
BEDFILE='../../../data/raw_data/meuleman_dhs_index/peaks/DHS_Index_and_Vocabulary.embryonic.full.narrowPeak.gz'
SOURCE='Homo_sapiens'
TARGET='Mus_musculus'
OUTDIR='../../../data/raw_data/meuleman_dhs_index/halper'
NAME=''

