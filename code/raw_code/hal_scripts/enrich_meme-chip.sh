#!/bin/bash
#SBATCH --partition=pfen1
#SBATCH --job-name=halliftover
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --error=logs/halliftover_%A.out.txt
#SBATCH --output=logs/halliftover_%A.out.txt

# load necessary packages
source ~/miniconda3/etc/profile.d/conda.sh
source ~/.bashrc
conda activate hal

function usage()
{
    echo "motif_enrichment.sh takes in a narrowpeakfile in one genome coordinate"
    echo "and uses meme-chip to run motif enrichment analyses."
    echo "Summits of narrowpeak files will be used to anchor peaks."
    echo "Example run call:"
    echo ""
    echo "sbatch motif_enrichment.sh -s sourceSpecies -t targetSpecies -o /dir/to/output/ -p myPeaks.narrowPeak "
    echo ""
    echo ""
    echo "Required parameters:"   
    echo "-p myPeaks.narrowPeak         narrowPeak file"
    echo "-o /dir/to/output/            output directory"
    echo "-s Mouse                      sourceSpecies in the Cactus file e.g. {Human, Mouse, Rhesus, etc}"
    echo "-t Human                      targetSpecies in the Cactus file e.g. {Human, Mouse, Rhesus, etc}"
    echo ""
    echo ""
    echo "Optional parameters:"   
    echo "-n myPeaks                    Mapped file name prefix after halliftover and Halper"
    echo "-c msaFile.hal                Cactus multi-species alignment file."
    echo "-f                            Whether to overwrite intermediate halliftover files."
    echo ""
}

# default cactus file
CACTUSFILE=/projects/pfenninggroup/machineLearningForComputationalBiology/alignCactus/10plusway-master.hal
OVERWRITE='FALSE'
NAME=''

# read in command line arguments
while [[ $1 != "" ]]; do
    case $1 in
        # the required parameters
        -p | --input-narrowPeak-file ) shift
                                NARROWPEAK=$1
                                echo "NarrowPeak file is $NARROWPEAK."
                                ;;
        -s | --source-species )      shift
                                SOURCE=$1
                                echo "Source species $SOURCE."
                                ;;
        -t | --target-species )      shift
                                TARGET=$1
                                echo "Target species $TARGET."
                                ;;
        -o | --output-dir )     shift
                                OUTDIR=$1
                                echo "Output dir is $OUTDIR."
                                ;;
        # some optional parameters
        -n | --name-base )      shift
                                NAME=$1
                                echo "Run name is $NAME."
                                ;;
        -c | --cactus-file )    shift
                                CACTUSFILE=$1
                                echo "Multi-species cactus alignment file is $CACTUSFILE."
                                ;;
        -f | --force-overwrite ) shift
                                OVERWRITE='TRUE'
                                echo "Overwrite intermediate halLiftover files."
                                ;;
        -h | --help )           usage
                                exit 1
                                ;;
        *)
        usage
        exit 1
    esac
    shift
done

########################################################################
# make output directory and working directory in scratch if do not exist
TMPDIR=/scratch/$USER/halLiftover
mkdir -p ${OUTDIR}; mkdir -p ${TMPDIR}

################################################################
# use the base name from the narrowPeak if name is not provided
if [[ ${NAME} == '' ]]
    then
    NAME=$(basename $NARROWPEAK | sed 's/.narrowPeak.gz//g;s/.narrowPeak//g')
fi
    

###################################
# convert narrowpeak to bed format
if [[ $(ls ${NARROWPEAK}| wc -l) == 0 ]]
    then
    echo "Input mouse narrowPeak ${NARROWPEAK} does not exist."
    exit 1
# check to see if narrowpeak exists
elif [[ ${NARROWPEAK} != '' ]]
    then
    # change narrowpeak to bed, give unique name, make sure has unique name
    echo "Converting narrowpeak file to bed file. Make sure has unique peak name in 4th column."
    BEDFILE=${TMPDIR}/${NAME}.tmp.bed
    if file --mime-type "${NARROWPEAK}" | grep -q gzip$;
        then zcat $NARROWPEAK | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $4, $5, "."}' > $BEDFILE
        else awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $4, $5, "."}' $NARROWPEAK > $BEDFILE
    fi
fi

##################################################################
# get summits from narrowpeak file or max cactus depth from cactus
SUMMITFILE=${TMPDIR}/${NAME}.summit.tmp.bed
echo "Extracting narrowPeak summits to halLiftover."
if file --mime-type "${NARROWPEAK}" | grep -q gzip$;
    # using the narrowpeak unique name
    then zcat $NARROWPEAK | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2+$10, $2+$10+1, $4, $5, "."}' > $SUMMITFILE
    else awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2+$10, $2+$10+1, $4, $5, "."}' $NARROWPEAK > $SUMMITFILE
fi

#####################################################
echo "Lifting ${SOURCE} summits to ${TARGET} genome."
HALLIFTEDSFILE=${OUTDIR}/${NAME}_${SOURCE}To${TARGET}_halLiftover.sFile.bed
if [[ ! -f ${HALLIFTEDSFILE}.gz || $OVERWRITE == 'TRUE' ]]
    then
    halLiftover $CACTUSFILE $SOURCE $SUMMITFILE $TARGET $HALLIFTEDSFILE
else
    echo "The file ${HALLIFTEDSFILE}.gz exist without permission to overwrite."
    echo "Unzipping and using the previous halliftover sFile."
    gunzip ${HALLIFTEDSFILE}.gz
fi

###################################################
echo "Lifting ${SOURCE} peaks to ${TARGET} genome."
HALLIFTEDTFILE=${OUTDIR}/${NAME}_${SOURCE}To${TARGET}_halLiftover.tFile.bed
if [[ ! -f ${HALLIFTEDTFILE}.gz || $OVERWRITE == 'TRUE' ]]
    then
    halLiftover $CACTUSFILE $SOURCE $BEDFILE $TARGET $HALLIFTEDTFILE
else
    echo "The file ${HALLIFTEDTFILE} exist without permission to overwrite."
    echo "Unzipping and using the previous halliftover tFile."
    gunzip ${HALLIFTEDTFILE}.gz
fi

########################################################################
echo "Apply HALPER to identify 1-1 orthologous in the ${TARGET} genome."
OUTFILE=${TMPDIR}/${NAME}_${SOURCE}To${TARGET}_orthologs.tmp.narrowPeak
python /home/bnphan/src/halLiftover-postprocessing/orthologFind.py \
-max_frac 1.5 -min_len 50 -protect_dist 5 -narrowPeak \
-tFile $HALLIFTEDTFILE \
-sFile $HALLIFTEDSFILE \
-qFile $BEDFILE \
-oFile $OUTFILE

#########################################
# sort hal-liftOvered mouse peaks to hg38
OUTFILE2=${OUTDIR}/${NAME}_${SOURCE}To${TARGET}_HALPER.narrowPeak
sort -k 1,1 -k2,2n $OUTFILE | uniq -u > $OUTFILE2
gzip --force $HALLIFTEDTFILE $HALLIFTEDSFILE $OUTFILE2
rm $BEDFILE $SUMMITFILE $OUTFILE

echo "Done."