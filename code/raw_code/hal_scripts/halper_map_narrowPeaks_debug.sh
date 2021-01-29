function cleanup()
{
    rm ${TMPDIR}/${NAME}*tmp*
}


# default cactus file
# CACTUSFILE=/projects/pfenninggroup/machineLearningForComputationalBiology/alignCactus/10plusway-master.hal
CACTUSFILE=/data/pfenninggroup/machineLearningForComputationalBiology/alignCactus/mam241/241-mammalian-2020v2.hal
OVERWRITE='FALSE'

cd /projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/code/raw_code/rheMac10_Stauffer_caudate

# read in command line arguments
BEDFILE='../../../data/raw_data/rheMac10/Stauffer_caudate/peak/Stauffer_caudate.Astro.GenBankRheMac8.narrowPeak.gz'
SOURCE='Macaca_mulatta'
TARGET='Homo_sapiens'
OUTDIR='../../../raw_data/rheMac10/Stauffer_caudate/peaks'
NAME=''

