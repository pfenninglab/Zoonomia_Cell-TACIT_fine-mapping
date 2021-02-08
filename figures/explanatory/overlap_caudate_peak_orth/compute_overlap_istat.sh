#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --time=1-0
#SBATCH --job-name=istat
#SBATCH --mem=12G
#SBATCH --error=logs/compute_overlap_istat_%A_%a.txt
#SBATCH --output=logs/compute_overlap_istat_%A_%a.txt

SCALE=100000

function usage()
{
    echo "compute_overlap_istat.sh takes in two bed-like files"
    echo "Both must have chr,start,end and can in either gzipped, or flat files."
    echo "Computes the overlap of the query bed file on the reference file using ISTAT."
    echo "Example run call:"
    echo ""
    echo "sbatch compute_overlap_istat.sh -q myQuery.narrowPeak -r myRef.bed.gz -o /path/to/out/dir"
    echo ""
    echo ""
    echo "Required parameters:"   
    echo "--query-bed-file  	-q FILENAME     bed or narrowPeak file, can be gzipped, {.bed/.bed.gz/.narrowpeak/.narrowpeak.gz}"
    echo "--reference-bed-file  -r FILENAME     bed or narrowPeak file, can be gzipped, {.bed/.bed.gz/.narrowpeak/.narrowpeak.gz}"
    echo "--output-dir      	-o /OUT/DIR/    path to main output directory"
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
        -q | --query-bed-file ) shift
                                QUERY=$1
                                echo "Query bed file is $QUERY."
                                ;;
        -r | --reference-bed-file )      shift
                                REF=$1
                                echo "Reference bed file is $REF."
                                ;;
        -o | --output-dir )     shift
                                OUTDIR=$1
                                echo "Output dir is $OUTDIR."
                                ;;
        *)
        usage
        exit 1
    esac
    shift
done

ISTAT=/home/bnphan/src/ISTAT/bin/istat
CHRSIZE=/home/bnphan/src/atac-seq-pipeline/genome/hg38/hg38.chrom.sizes
mkdir -p $OUTDIR

## create scratch space to work in ##
TMP_LABEL="tmp$(date +%s)"
TMPDIR=/scratch/bnphan/istat/${TMP_LABEL}
mkdir -p ${TMPDIR}; cd $TMPDIR

## copy query bed file
QUERYNAME=$(basename ${QUERY} .narrowPeak.gz)
QUERYBED=${TMPDIR}/${QUERYNAME}.query.bed
if file --mime-type "${QUERY}" | grep -q gzip$ ;
    then zcat $QUERY | sort -k 1,1 -k2,2n | cut -f1-3 > $QUERYBED
    else cat $QUERY | sort -k 1,1 -k2,2n | cut -f1-3 > $QUERYBED
fi

## copy reference bed file
REFNAME=$(basename ${REF} .narrowPeak.gz)
REFBED=${TMPDIR}/${REFNAME}.reference.bed
if file --mime-type "${REF}" | grep -q gzip$ ;
    then zcat $REF | sort -k 1,1 -k2,2n | cut -f1-3 > $REFBED
    else cat $REF | sort -k 1,1 -k2,2n | cut -f1-3 > $REFBED
fi

OUTPUT_FILE=${TMPDIR}/${QUERYNAME}_#_${REFNAME}.istat.txt

### Run iSTAT ###
${ISTAT} ${QUERYBED} ${REFBED} ${CHRSIZE} ${SCALE} ${OUTPUT_FILE} s
rsync --remove-source-files -Paq ${TMPDIR}/${QUERYNAME}.${REFNAME}.istat.txt ${OUTDIR}
rm -r ${TMPDIR}