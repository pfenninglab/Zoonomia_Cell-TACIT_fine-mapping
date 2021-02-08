import pandas as pd
import pybedtools
import os
import numpy as np

GWASDIR='/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments'
POLYFUNDIR=f'{GWASDIR}/polyfun/baselineLF2.2.UKB'
SETDIR='/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate'
BEDDIR=f'{SETDIR}/data/raw_data/snp_annotations/bed'
OUTDIR=f'{SETDIR}/data/raw_data/snp_annotations/halper_zoo'
ZOONOMIADIR=f'{SETDIR}/data/tidy_data/Zoonomia_data'

os.chdir(f'{SETDIR}/code/raw_code/snp_annotations')
os.system(f'mkdir -p {BEDDIR} {OUTDIR} logs')
# os.getcwd()

zoo_df = pd.read_csv(f'{ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv', sep = '\t')
# zoo_df = zoo_df.sort_values(by='Zoonomia Index')
zoo_df = zoo_df[zoo_df['Species'] != 'Homo_sapiens']
zoo_df = zoo_df.sort_index(ascending=False)
TARGETS=zoo_df['Species'].str.cat(sep=',')
SOURCE='Homo_sapiens'
SLURM_CALL=f'sbatch --mem 4G -p pfen1 -w compute-1-39 {SETDIR}/code/raw_code/hal_scripts/halper_map_peak_orthologs.sh'

# prepare empty dataframe to stash SNPs w/o db150 or db151 snps
notFound = pd.DataFrame()

for CHR in range(1,23):
	ANNOTFN=f"{BEDDIR}/UKB_imputed_snps.GRCh38.{CHR}.bed.gz"
	BEDFILE150=f'{GWASDIR}/1000G_ALL_Phase3_hg38_files/vcf/b150_GRCh38p13_snps_bed/bed_chr_{CHR}.bed.gz'
	BEDFILE151=f'{GWASDIR}/1000G_ALL_Phase3_hg38_files/vcf/b151_GRCh38p7_snps_bed/bed_chr_{CHR}.bed.gz'
	if not os.path.isfile(ANNOTFN):
		print(f'Extracting SNPs from UKB chr{CHR}.')
		df_annot = pd.read_parquet(f'{POLYFUNDIR}/weights.UKB.{CHR}.l2.ldscore.parquet')
		df_annot['NAME'] = df_annot['SNP'].astype(str) +'_' + df_annot['A1'].astype(str) +'_' +df_annot['A2'].astype(str)
		# read in the GRCh38 SNP coordinates
		df_b150 = pd.read_csv(BEDFILE150, names = ['SEQ','START','END','SNP'], index_col=False, 
		skiprows = 1, compression = 'gzip', sep = '\t')
		df_b151 = pd.read_csv(BEDFILE151, names = ['SEQ','START','END','SNP'], index_col=False, 
		skiprows = 1, compression = 'gzip', sep = '\t')
		df_dbsnp = pd.concat([df_b150, df_b151], axis=0).drop_duplicates(subset=['SNP', 'START','END'])
		# get the hg38 bed coordinates of UKB snps, drop missing in UKB
		df_out = pd.merge(df_annot, df_dbsnp, how = 'left', on='SNP')
		notFound = notFound.append(df_out[~df_out['SEQ'].notna()])
		df_out = df_out[df_out['SEQ'].notna()]
		df_out.START = df_out.START.astype(int)
		df_out.END = df_out.END.astype(int)
		indZeroWidth = df_out.START == df_out.END
		df_out.loc[indZeroWidth, ('END')] = df_out.loc[indZeroWidth, ('START')] + 1
		np.sum(df_out.START ==df_out.END)
		# read in the GRCh38 SNP coordinates, by rsID
		# df_out.to_csv(ANNOTFN, columns = ['SEQ','START','END','NAME'],
		# header = False, index=False, sep = '\t', compression="gzip")

for CHR in range(1,23):
	ANNOTFN=f"{BEDDIR}/UKB_imputed_snps.GRCh38.{CHR}.bed.gz"
	THECALL=f'{SLURM_CALL} -s {SOURCE} -t {TARGETS} -o {OUTDIR} -b {ANNOTFN} --snp'
	print(f'Halper UKB SNPs on {CHR} from {TARGETS}.')
	# os.system(THECALL)


ANNOTFN2=f"{BEDDIR}/UKB_imputed_snps.hg19.unmapped.bed"
ANNOTFN3=f"{BEDDIR}/UKB_imputed_snps.GRCh38.unmappedlifted.bed"
if not os.path.isfile(f'{ANNOTFN3}.gz'):
	df_toLiftOver = notFound.copy()
	df_toLiftOver['START'] = df_toLiftOver['BP'] - 1
	df_toLiftOver['SEQ'] = 'chr' + df_toLiftOver['CHR'].astype('str').str.strip('chr')
	df_toLiftOver.to_csv(ANNOTFN2, columns = ['SEQ','START','BP','NAME'],
			header = False, index=False, sep = '\t')
	CHAIN='/home/bnphan/resources/liftOver_chainz/hg19ToHg38.over.chain'
	THECALL1=f'liftOver {ANNOTFN2} {CHAIN} {ANNOTFN3} {ANNOTFN3}.unmapped'
	THECALL2=f'gzip {ANNOTFN3}'
	os.system(THECALL1)
	os.system(THECALL2)

THECALL3=f'{SLURM_CALL} -s {SOURCE} -t {TARGETS} -o {OUTDIR} -b {ANNOTFN3}.gz --snp'
os.system(THECALL3)




