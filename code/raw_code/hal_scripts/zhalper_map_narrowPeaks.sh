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
    echo "halper_map_peak_orthologs.sh takes in a bed/narrowpeak file in one genome coordinate"
    echo "and uses halLiftover + HALPER to map to the orthologous regions in a target genome."
    echo "Summits of narrowpeak files will be used to anchor peaks."
    echo "Example run call:"
    echo ""
    echo "sbatch halper_map_peak_orthologs.sh -b myPeaks.narrowPeak -o /dir/to/output/ -s sourceSpecies -t targetSpecies "
    echo ""
    echo ""
    echo "Required parameters:"   
    echo "-b FILENAME   bed or narrowPeak file, can be gzipped, {.bed/.bed.gz/.narrowpeak/.narrowpeak.gz}"
    echo "-o /OUT/DIR/  path to main output directory"
    echo "-s SPECIES    sourceSpecies in the Cactus file e.g. {Human, Homo_sapiens, Rhesus, etc}"
    echo "-t SPECIES    targetSpecies in the Cactus file e.g. {Human, Mouse, Rhesus, etc}"
    echo ""
    echo ""
    echo "Optional parameters:"   
    echo "-n NAME       Mapped file name prefix after halliftover and Halper"
    echo "-c HALFILE    Cactus multi-species alignment file {msaFile.hal}"
    echo "-f            Whether to overwrite intermediate halliftover files."
    echo ""
    echo "Note: run `halStats file.hal` to get list of genomes species."
    echo ""
}

function cleanup()
{
    rm ${TMPDIR}/${NAME}*tmp*
}

# default cactus file
# CACTUSFILE=/projects/pfenninggroup/machineLearningForComputationalBiology/alignCactus/10plusway-master.hal
CACTUSFILE=/data/pfenninggroup/machineLearningForComputationalBiology/alignCactus/mam241/241-mammalian-2020v2.hal
OVERWRITE='FALSE'
NAME='';

# read in command line arguments
while [[ $1 != "" ]]; do
    case $1 in
        # the required parameters
        -b | --input-bed-file ) shift
                                BEDFILE=$1
                                echo "Bed file is $BEDFILE."
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

############################################################
# check required parameters, and some sensibilites. ########
if [[ -z "$SOURCE" ]];
    then echo 'Please set source species e.g. {Human, Mouse, Rhesus, etc}.'; exit 1
elif [[ -z "$TARGET" ]];
    then echo 'Please set target species e.g. {Human, Mouse, Rhesus, etc}.'; exit 1
elif [[ -z "$OUTDIR" ]];
    then echo 'Please set output directory.'; exit 1
elif [[ "$SOURCE" == "$TARGET" ]];
    then echo "Source and target species cannot be the same."; exit 1
elif [[ $(ls ${BEDFILE}| wc -l) == 0 ]]
    then echo "Input bed/narrowpeak file does not exist:${BEDFILE}"; exit 1
elif [[ ${NAME} == '' ]]
    then NAME=$(basename $BEDFILE | sed 's/.gz$//g;s/.narrowPeak$//g;s/.bed$//g')
fi

########################################################################
# make output directory and working directory in scratch if do not exist
TMPDIR=/scratch/$USER/halLiftover
mkdir -p ${OUTDIR} ${TMPDIR}

# unzip gzipped bed file
INPUTBED=${TMPDIR}/${NAME}.input.tmp.bed
if file --mime-type "${BEDFILE}" | grep -q gzip$ ;
    then zcat $BEDFILE > $INPUTBED
    else cat $BEDFILE > $INPUTBED
fi

##########################################
UNIQUEBED=${TMPDIR}/${NAME}.unique.tmp.bed
# check if duplicate peaks
if [[ $(awk '++A[$1,$2,$3,$10] > 1 { print "true"; exit 1 }' $INPUTBED) ]]; then
    echo "Repeat chr:start-end:summits detected. Remove duplicate peaks before proceeding."; cleanup; exit 1
# check if basic 3 column bed file w/o name column
elif [[ $(awk '{print NF; exit}' ${INPUTBED}) -lt 4 ]]; 
    then echo "Bed file doesn't have name column. Adding"
    echo "USING CHR:START-END in the NAME column."
    awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $1":"$2-$3"," 0, "."}' $INPUTBED > $UNIQUEBED    
# name column found, check if duplicate names
elif [[ $(awk '++A[$4] > 1 { print "true"; exit 1 }' $INPUTBED) ]]; then
    echo "Non-unique bed peak names detected. Giving unique names now."
    echo "Bed file doesn't have name column. Adding"
    if [[ $(awk '{print NF; exit}' ${INPUTBED}) == 10 ]]; 
    # use summit if there's a 10th column (assume narrowpeak file)
        then echo "Appending CHR:START-END:SUMMIT to NAME column."
        awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $4":"$1":"$2-$3":"$10, $5, $6}' $INPUTBED > $UNIQUEBED
    else # if there isn't a narrowpeak summit column
        echo "Appending CHR:START-END to NAME column."
        awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $4":"$1":"$2-$3"," $5, $6}' $INPUTBED > $UNIQUEBED
    fi
else 
    echo "Bed peak names are unique; moving onwards."
    cat $INPUTBED > $UNIQUEBED
fi

#################################################################
# make sure the scores column is numeric, strand column is +|-|.
echo "Formatting bed score and strand columns."
awk 'BEGIN {FS="\t"; OFS="\t"} {! $5 ~ /^[[:digit:]]+$/} {$5=0} {print;}' $UNIQUEBED > ${TMPDIR}/${NAME}.unique2.tmp.bed
awk 'BEGIN {FS="\t"; OFS="\t"} {! $6 ~ /+|-|./} {$6="."} {print;}' ${TMPDIR}/${NAME}.unique2.tmp.bed > ${TMPDIR}/${NAME}.unique3.tmp.bed
mv ${TMPDIR}/${NAME}.unique3.tmp.bed $UNIQUEBED; rm ${TMPDIR}/${NAME}.unique2.tmp.bed

##################################################################
# get summits from narrowpeak file or max cactus depth from cactus
SUMMITFILE=${TMPDIR}/${NAME}.summit.tmp.bed
echo "Extracting bed/narrowpeak summits to halLiftover."
if [[ $(awk '{print NF; exit}' ${UNIQUEBED}) -lt 3 ]]; 
    then echo "Too few columns to be a bed file."; cleanup; exit 1
elif [[ $(awk '{print NF; exit}' ${UNIQUEBED}) == 10 ]]; 
    # 10 columns, assume this is a narrowpeak
    then echo "Summits detected. Using the summits."
    awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2+$10, $2+$10+1, $4, $5, $6}' $UNIQUEBED > $SUMMITFILE
else [[ $(awk '{print NF; exit}' ${UNIQUEBED}) -gt 3 ]]
    # take the mean value between start and end of each peak
    echo "No summits found, taking the mean between start and end as the summits."
    awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, int(($2+$3)/2), int(($2+$3)/2)+1, $4, $5, $6}' $UNIQUEBED > $SUMMITFILE
fi

######################################################
echo "Lifting ${SOURCE} summits to ${TARGET} genome."
HALLIFTEDSFILE=${OUTDIR}/${NAME}.${SOURCE}To${TARGET}.halLiftover.sFile.bed
if [[ ! -f ${HALLIFTEDSFILE}.gz || $OVERWRITE == 'TRUE' ]]
    then halLiftover $CACTUSFILE $SOURCE $SUMMITFILE $TARGET $HALLIFTEDSFILE
else
    echo "The file ${HALLIFTEDSFILE}.gz exist without permission to overwrite."
    echo "Unzipping and using the previous halliftover sFile."
    gunzip ${HALLIFTEDSFILE}.gz
fi

###################################################
echo "Lifting ${SOURCE} peaks to ${TARGET} genome."
HALLIFTEDTFILE=${OUTDIR}/${NAME}.${SOURCE}To${TARGET}.halLiftover.tFile.bed
PEAKFILE=${TMPDIR}/${NAME}.peak.tmp.bed
cut -f 1-6 $UNIQUEBED > $PEAKFILE # 6 column bed format
if [[ ! -f ${HALLIFTEDTFILE}.gz || $OVERWRITE == 'TRUE' ]]
    then halLiftover $CACTUSFILE $SOURCE $PEAKFILE $TARGET $HALLIFTEDTFILE
else
    echo "The file ${HALLIFTEDTFILE}.gz exists without permission to overwrite."
    echo "Unzipping and using the previous halliftover tFile."
    gunzip ${HALLIFTEDTFILE}.gz
fi

########################################################################
echo "Apply HALPER to identify 1-1 orthologous in the ${TARGET} genome."
OUTFILE=${TMPDIR}/${NAME}.${SOURCE}To${TARGET}.orthologs.tmp.narrowPeak
python /home/bnphan/src/halLiftover-postprocessing/orthologFind.py \
-max_frac 1.5 -min_len 50 -protect_dist 10 -narrowPeak \
-tFile $HALLIFTEDTFILE -sFile $HALLIFTEDSFILE \
-qFile $UNIQUEBED -oFile $OUTFILE

echo "Mapped $(wc -l $OUTFILE | cut -d ' ' -f1) peaks out of $(wc -l $INPUTBED | cut -d ' ' -f1) total peaks."

#########################################
# sort hal-liftOvered mouse peaks to hg38
OUTFILE2=${OUTDIR}/${NAME}.${SOURCE}To${TARGET}.HALPER.narrowPeak
sort -k 1,1 -k2,2n $OUTFILE | uniq -u > $OUTFILE2
gzip --force $HALLIFTEDTFILE $HALLIFTEDSFILE $OUTFILE2
echo "Done."; cleanup


