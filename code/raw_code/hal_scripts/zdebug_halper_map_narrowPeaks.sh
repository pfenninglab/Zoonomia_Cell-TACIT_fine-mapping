# read in command line arguments
cd /projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/code/raw_code/snp_annotations

BEDFILE='../../../data/raw_data/snp_annotations/bed/UKB_imputed_snps.GRCh38.22.bed.gz'
SOURCE='Homo_sapiens'
TARGETS='Macaca_mulatta,Mus_musculus'
OUTDIR='../../../data/raw_data/snp_annotations/halper_zoo'
NAME=''
MIN_LEN=0
PROTECT_DIST=0
SNP='TRUE'
TARGET='Macaca_mulatta'


# read in command line arguments
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/mm10_BICCN_mouse_caudoputamen
DATADIR=${SETDIR}/data/raw_data/mm10/Corces_2020
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments

cd $CODEDIR
BEDFILE='/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/raw_data/mm10/Mouse_cSNAIL_D1D2/peak/Pfenning_Cpu.MSN_D2.narrowPeak.gz'
SOURCE='Mus_musculus'
TARGETS='Macaca_mulatta,Homo_sapiens'
TARGET='Macaca_mulatta'
OUTDIR='/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/raw_data/mm10/Mouse_cSNAIL_D1D2/halper'
NAME=''


QUERY=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/raw_data/hg38/Corces_2020/peak/Corces2020_caudate_hgMmOrth.MSN_D1.narrowPeak.gz
REF=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/raw_data/rheMac10/Stauffer_caudate/halper/Stauffer_caudate.MSN_D1.GenBankRheMac8.Macaca_mulattaToHomo_sapiens.HALPER.narrowPeak.gz


