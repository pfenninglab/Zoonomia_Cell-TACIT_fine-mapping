# directories 
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/snp_annotations
DATADIR=${SETDIR}/data/raw_data/snp_annotations
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
SNPDIR=${SETDIR}/data/tidy_data/ldsc_gwas/1000G_ALL_Phase3_hg38_files/vcf
cd $CODEDIR; mkdir -p $SNPDIR


##################
# make plink bim files from HapMap3 SNPs sans HLA to bed format
for BIM in ${GWASDIR}/1000G_ALL_Phase3_hg38_files/plink/*HM3*.bim; do
SNPPOS=${SNPDIR}/$(basename $BIM .bim).snp.bed.gz
# .bim in 1-base, bed in 0-base
# adjust for multi-nucleotide polymoriphisms, sizing w/ A1 allele
awk -v OFS='\t' '{print "chr"$1, $4-1, $4 + length($5) - 1, $2}' $BIM | gzip > $SNPPOS
done

cat ${SNPDIR}/*AFR*.snp.bed.gz | sort -k1,1 -k2,2n | gzip > ${SNPDIR}/1000G.AFR.HM3.snp.bed.gz
cat ${SNPDIR}/*EUR*.snp.bed.gz | sort -k1,1 -k2,2n | gzip > ${SNPDIR}/1000G.EUR.HM3.snp.bed.gz
rm ${SNPDIR}/*.snp.bed.gz

#############################################
# halper the HM3 SNPs across Zoonomia genomes 
SOURCE='Homo_sapiens'
OUTDIR=${DATADIR}/haper_snps

for IND in {1..243}; do
for BEDFILE in ${SNPDIR}/*.snp.bed.gz; do
TARGET=$(awk -F'\t' -v IND=$IND 'FNR == (IND + 1 ) {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv)
NAME=$(basename $BEDFILE | sed 's/.gz$//g;s/.narrowPeak$//g;s/.bed$//g')
if [[ $TARGET != 'Homo_sapiens' && ! -f ${OUTDIR}/${NAME}.${SOURCE}To${TARGET}.HALPER.narrowPeak.gz ]]; then
sbatch -p pfen1 -w compute-1-40 \
${CODEDIR}/../hal_scripts/halper_map_peak_orthologs.sh \
-s ${SOURCE} -t ${TARGET} -o ${OUTDIR} -b ${BEDFILE}
fi
done
done


for CHR in {1..22}; do
zcat ${GWASDIR}/1000G_ALL_Phase3_hg38_files/vcf/All_20180418_b151_GRCh38p7.vcf.gz | \
sed -e 's/chr//' | awk -v CHR=$CHR '{OFS="\t"; if (!/^#/ && CHR == $1) {print "chr"$1,$2-1,$2 + length($4) - 1, $3"_"$4"_"$5}}' | \
gzip > ${SNPDIR}/All_20180418_b151_GRCh38p7.${CHR}.bed.gz
done


for CHR in {1..22}; do
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/BED/bed_chr_${CHR}.bed.gz \
-P ${GWASDIR}/1000G_ALL_Phase3_hg38_files/vcf/b151_GRCh38p7_snps_bed

wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/BED/bed_chr_${CHR}.bed.gz \
-P ${GWASDIR}/1000G_ALL_Phase3_hg38_files/vcf/b150_GRCh38p13_snps_bed
done

