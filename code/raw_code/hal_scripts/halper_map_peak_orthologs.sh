#!/bin/bash
#SBATCH --partition=pfen1
#SBATCH --job-name=halliftover
#SBATCH --ntasks-per-core=1
#SBATCH --mem=10G
#SBATCH --error=logs/halliftover_%A.out.txt
#SBATCH --output=logs/halliftover_%A.out.txt

# load necessary packages
source ~/miniconda3/etc/profile.d/conda.sh
source ~/.bashrc

# CACTUSFILE=/projects/pfenninggroup/machineLearningForComputationalBiology/alignCactus/10plusway-master.hal
conda activate hal
CACTUSFILE=/data/pfenninggroup/machineLearningForComputationalBiology/alignCactus/mam241/241-mammalian-2020v2.hal
OVERWRITE='FALSE'; NAME=''; SNP=''; MIN_LEN=50; PROTECT_DIST=10 # for halper param
TMP_LABEL="tmp$(date +%s)"

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
    echo "--input-bed-file  -b FILENAME     bed or narrowPeak file, can be gzipped, {.bed/.bed.gz/.narrowpeak/.narrowpeak.gz}"
    echo "--output-dir      -o /OUT/DIR/    path to main output directory"
    echo "--source-species  -s SPECIES      sourceSpecies in the Cactus file e.g. {Homo_sapiens, Mus_musculus, Macaca_mulatta}"
    echo "--target-species  -t SPECIES1,SPECIES2"
    echo "                                  comma separated list of targetSpecies in the Cactus file e.g. {Human, Mouse, Rhesus, etc}"
    echo ""
    echo ""
    echo "Optional parameters:"   
    echo "--name-base       -n NAME         Mapped file name prefix after halliftover and Halper"
    echo "--cactus-file     -c HALFILE      Cactus multi-species alignment file {msaFile.hal}"
    echo "--force-overwrite -f              Whether to overwrite intermediate halliftover files."
    echo "--snp                             Map SNPs rather than peaks (change min peak length=1/protect dist=0)."
    echo ""
    echo ""
    echo "Note: genomes are whatever you get from running:"
    echo "  halStats --genomes file.hal"
    echo ""
}

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
                                TARGETS=$1
                                echo "Target species $TARGETS."
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
        --snp )                 shift
                                SNP='TRUE';
                                echo "Mapping SNPs rather than peaks. Setting min length=1."
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

function check_params()
{
    ############################################################
    # check required parameters, and some sensibilites. ########
    if [[ -z "$SOURCE" ]];
        then echo 'Please set source species e.g. {Human, Mouse, Rhesus, etc}.'; exit 1
    elif [[ -z "$TARGETS" ]];
        then echo 'Please set target species e.g. {Human, Mouse, Rhesus, etc}.'; exit 1
    elif [[ -z "$OUTDIR" ]];
        then echo 'Please set output directory.'; exit 1
    elif [[ "$SOURCE" == "$TARGETS" ]];
        then echo "Source and target species cannot be the same."; exit 1
    elif [[ $(ls ${BEDFILE}| wc -l) == 0 ]]
        then echo "Input bed/narrowpeak file does not exist:${BEDFILE}"; exit 1
    elif [[ ${NAME} == '' ]]
        then NAME=$(basename $BEDFILE | sed 's/.gz$//g;s/.narrowPeak$//g;s/.bed$//g')
    fi
    TARGETS=$(echo $TARGETS | tr "," "\n")
}

function check_bed()
{
    # unzip gzipped bed file
    INPUTBED=${TMPDIR}/${NAME}.input.${SOURCE}.${TMP_LABEL}.bed
    if file --mime-type "${BEDFILE}" | grep -q gzip$ ;
        then zcat $BEDFILE > $INPUTBED
        else cat $BEDFILE > $INPUTBED
    fi
}

function format_bed()
{
    ##########################################
    UNIQUEBED=${TMPDIR}/${NAME}.unique.${SOURCE}.${TMP_LABEL}.bed
    # check if basic 3 column bed file w/o name column
    if [[ $(awk '{print NF; exit}' ${INPUTBED}) -lt 4 ]]; 
        then echo "Bed file doesn't have name column. Adding"
        echo "USING CHR:START-END in the NAME column."
        awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $1 ":" $2 "-" $3 ":" int(($3 - $2) / 2), "0", "."}' $INPUTBED > $UNIQUEBED    
    # name column found, check if duplicate names
    elif [[ $(awk '++A[$4] > 1 { print "true"; exit 1 }' $INPUTBED) ]]; then
        echo "Non-unique bed peak names detected. Giving unique names now."
        echo "Bed file doesn't have name column. Adding"
        if [[ $(awk '{print NF; exit}' ${INPUTBED}) == 10 ]]; 
        # use summit if there's a 10th column (assume narrowpeak file)
            then echo "Appending CHR:START-END:SUMMIT to NAME column."
            awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $1 ":" $2 "-" $3 ":" $10, $5, $6}' $INPUTBED > $UNIQUEBED
        else # if there isn't a narrowpeak summit column
            echo "Appending CHR:START-END to NAME column."
            awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $1 ":" $2 "-" $3, $5, $6}' $INPUTBED > $UNIQUEBED
        fi
    else 
        echo "Bed peak names are unique; moving onwards."
        cat $INPUTBED > $UNIQUEBED
    fi

    #################################################################
    # make sure the scores column is numeric, strand column is +|-|.
    echo "Formatting bed score and strand columns."
    awk 'BEGIN {FS="\t"; OFS="\t"} {! $5 ~ /^[[:digit:]]+$/} {$5=0} {print;}' $UNIQUEBED \
        > ${TMPDIR}/${NAME}.unique2.${SOURCE}.${TMP_LABEL}.bed
    awk 'BEGIN {FS="\t"; OFS="\t"} {! $6 ~ /+|-|./} {$6="."} {print;}' \
        ${TMPDIR}/${NAME}.unique2.${SOURCE}.${TMP_LABEL}.bed > ${TMPDIR}/${NAME}.unique3.${SOURCE}.${TMP_LABEL}.bed
    mv ${TMPDIR}/${NAME}.unique3.${SOURCE}.${TMP_LABEL}.bed $UNIQUEBED; rm ${TMPDIR}/${NAME}.unique2.${SOURCE}.${TMP_LABEL}.bed

    ######################################
    # get the 6-column bed for broak peak
    SIMPLEBED=${TMPDIR}/${NAME}.6col.${SOURCE}.${TMP_LABEL}.bed; cut -f1-6 $UNIQUEBED > $SIMPLEBED
}

function get_summits()
{
    ##################################################################
    # get summits from narrowpeak file or max cactus depth from cactus
    SUMMITFILE=${TMPDIR}/${NAME}.summit.${SOURCE}.${TMP_LABEL}.bed
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
}

function get_cactus()
{
    ########################################################################
    # copy over cactus file to speed up the IO during hal-liftover, need > 1Tb
    TMPDIR=/scratch/$USER/halLiftover; mkdir -p ${OUTDIR} ${TMPDIR}
    rsync -Paq $CACTUSFILE $TMPDIR
}

function lift_summits()
{
    ######################################################
    echo "Lifting ${SOURCE} summits to ${TARGET} genome."
    HALLIFTEDSFILE=${NAME}.${SOURCE}To${TARGET}.halLiftover.sFile.bed
    if [[ ! -f ${OUTDIR}/${HALLIFTEDSFILE}.gz || $OVERWRITE == 'TRUE' ]]
        # hal-liftover the summits de novo
        then halLiftover $TMPDIR/$(basename $CACTUSFILE) $SOURCE $SUMMITFILE $TARGET ${TMPDIR}/$HALLIFTEDSFILE
    else
        # using previously hal-liftover summits
        echo "The file ${HALLIFTEDSFILE}.gz exist without permission to overwrite."
        echo "Unzipping and using the previous halliftover sFile."
        gunzip -c ${OUTDIR}/${HALLIFTEDSFILE}.gz > ${TMPDIR}/${HALLIFTEDSFILE}
    fi
}

function lift_peaks()
{
    ###################################################
    echo "Lifting ${SOURCE} peaks to ${TARGET} genome."
    HALLIFTEDTFILE=${NAME}.${SOURCE}To${TARGET}.halLiftover.tFile.bed
    if [[ ! -f ${OUTDIR}/${HALLIFTEDTFILE}.gz || $OVERWRITE == 'TRUE' ]]
        # hal-liftover the peaks de novo
        then halLiftover $TMPDIR/$(basename $CACTUSFILE) $SOURCE $SIMPLEBED $TARGET ${TMPDIR}/$HALLIFTEDTFILE
    else
        # using previously hal-liftover peaks
        echo "The file ${HALLIFTEDTFILE}.gz exists without permission to overwrite."
        echo "Unzipping and using the previous halliftover tFile."
        gunzip -c ${OUTDIR}/${HALLIFTEDTFILE}.gz > ${TMPDIR}/${HALLIFTEDTFILE}
    fi
}

function run_halper()
{    
    ########################################################################
    echo "Apply HALPER to identify 1-1 orthologous in the ${TARGET} genome."
    OUTFILE=${TMPDIR}/${NAME}.${SOURCE}To${TARGET}.orthologs.${SOURCE}To${TARGET}.${TMP_LABEL}.narrowPeak
    python /home/bnphan/src/halLiftover-postprocessing/orthologFind.py \
    -max_frac 1.5 -min_len $MIN_LEN -protect_dist $PROTECT_DIST \
    -tFile ${TMPDIR}/$HALLIFTEDTFILE -sFile ${TMPDIR}/$HALLIFTEDSFILE \
    -narrowPeak -qFile $UNIQUEBED -oFile $OUTFILE
    echo "Mapped $(wc -l $OUTFILE | cut -d ' ' -f1) peaks out of $(wc -l $INPUTBED | cut -d ' ' -f1) total peaks."

    #################################################################
    # sort hal-liftOvered mouse peaks to hg38 and copy to destination
    # if [[ ${SNP} == "TRUE" ]]; then rm ${TMPDIR}/$HALLIFTEDTFILE; fi
    OUTFILE2=${TMPDIR}/${NAME}.${SOURCE}To${TARGET}.HALPER.narrowPeak
    sort -k 1,1 -k2,2n $OUTFILE | uniq -u > $OUTFILE2
    gzip --force ${TMPDIR}/$HALLIFTEDTFILE ${TMPDIR}/$HALLIFTEDSFILE $OUTFILE2
    rsync --remove-source-files -Paq ${TMPDIR}/${NAME}.${SOURCE}To${TARGET}*.gz ${OUTDIR}/
    rm ${TMPDIR}/${NAME}*${SOURCE}To${TARGET}.${TMP_LABEL}*
    echo "Done."; 
}

function map_snps()
{
    ####################################################################
    echo "Only the ${SOURCE} bed files of the SNPs to ${TARGET} genome."
    HALLIFTEDTFILE=${NAME}.${SOURCE}To${TARGET}.halLiftover.tFile.bed
    OUTFILE=${TMPDIR}/${NAME}.${SOURCE}To${TARGET}.snp.bed
    if [[ ! -f ${OUTDIR}/${HALLIFTEDTFILE}.gz || $OVERWRITE == 'TRUE' ]]
        # hal-liftover the peaks de novo
        then halLiftover $TMPDIR/$(basename $CACTUSFILE) $SOURCE $SIMPLEBED $TARGET $OUTFILE
        sort -k 1,1 -k2,2n $OUTFILE | uniq -u > ${TMPDIR}/$HALLIFTEDTFILE
        rm $OUTFILE
    else
        # using previously hal-liftover peaks
        echo "The file ${HALLIFTEDTFILE}.gz exists without permission to overwrite."
        echo "Unzipping and using the previous halliftover tFile."
        gunzip -c ${OUTDIR}/${HALLIFTEDTFILE}.gz > ${TMPDIR}/${HALLIFTEDTFILE}
    fi
    echo "Mapped $(wc -l ${TMPDIR}/$HALLIFTEDTFILE | cut -d ' ' -f1) out of $(wc -l $INPUTBED | cut -d ' ' -f1) total SNPs."

    ######################################
    gzip --force ${TMPDIR}/$HALLIFTEDTFILE
    rsync --remove-source-files -Paq ${TMPDIR}/$HALLIFTEDTFILE.gz ${OUTDIR}/
}

function cleanup()
{
    rm ${TMPDIR}/*${TMP_LABEL}*
}

##############################
# prepare bed files for input
check_params
check_bed
format_bed
get_cactus

##########################################
# perform liftover for each target species
for TARGET in $TARGETS; do
    if [[ "${SNP}" == 'TRUE' ]]; then map_snps
    else get_summits; lift_summits; lift_peaks; run_halper; fi
done

##################
# final clean up
cleanup
