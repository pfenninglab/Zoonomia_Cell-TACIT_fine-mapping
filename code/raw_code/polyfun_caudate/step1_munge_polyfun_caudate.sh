#!/bin/bash
SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate'
GWASDIR='/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas/unmunged_sumstats'
VCFDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_ALL_Phase3_hg19_files/vcf
DATADIR=${SETWD}/data/raw_data/polyfun_caudate
CODEDIR=${SETWD}/code/raw_code/polyfun_caudate
POLYFUNDIR='/home/bnphan/src/polyfun'

cd $CODEDIR; mkdir -p $DATADIR $DATADIR/munged;
sed -i 's/functional_polyfun_Zoonomia_conservation/polyfun_caudate/g' $CODEDIR/*

source activate polyfun
rsync -Pav /home/bnphan/projects/Addiction_MPRA_2021/data/raw_data/polyfun/munged $DATADIR

###################
## 1) Degen traits 
LABEL=AD-Kunkle_2019
SUMSTATS=${GWASDIR}/Alzheimers_GWAS/Kunkle_etal_Stage1_minProt.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--n 63926 --min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi

LABEL=AD-Jansen_2019
SUMSTATS=${GWASDIR}/Alzheimers_GWAS/AD_sumstats_Jansenetal.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


###################
## 2) Psych traits 
if [[ ! -f ${DATADIR}/munged/Cross_disorder-PGCCDG_2019.parquet ]]; then
python ${POLYFUNDIR}/munge_polyfun_sumstats.py --min-info 0.6 --min-maf 0.001 \
--sumstats ${GWASDIR}/PGC_ADHD_BP_OCD_MDD_SCZ3_PTSD/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt.gz \
--out ${DATADIR}/munged/Cross_disorder-PGCCDG_2019.parquet; fi


if [[ ! -f ${DATADIR}/munged/ADHD-Demontis_2018.parquet ]]; then
zcat ${GWASDIR}/PGC_ADHD_BP_OCD_MDD_SCZ3_PTSD/adhd_jul2017_sumstats.gz | grep -Ev $'NA' | \
awk '!seen[$2]++' | gzip > ${GWASDIR}/PGC_ADHD_BP_OCD_MDD_SCZ3_PTSD/adhd_jul2017_sumstats_nodup.gz 
python ${POLYFUNDIR}/munge_polyfun_sumstats.py --min-info 0.6 --min-maf 0.001 \
--n 55374 --sumstats ${GWASDIR}/PGC_ADHD_BP_OCD_MDD_SCZ3_PTSD/adhd_jul2017_sumstats_nodup.gz \
--out ${DATADIR}/munged/ADHD-Demontis_2018.parquet; fi


if [[ ! -f ${DATADIR}/munged/Anorexia-Watson_2019.parquet ]]; then
python ${POLYFUNDIR}/munge_polyfun_sumstats.py --min-info 0.6 --min-maf 0.001 \
--sumstats ${GWASDIR}/PGC_ADHD_BP_OCD_MDD_SCZ3_PTSD/pgcAN2.2019-07.trun.tsv.gz \
--out ${DATADIR}/munged/Anorexia-Watson_2019.parquet; fi


if [[ ! -f ${DATADIR}/munged/Bipolar-Mullins_2021.parquet ]]; then
python ${POLYFUNDIR}/munge_polyfun_sumstats.py --min-info 0.6 --min-maf 0.001 \
--sumstats ${GWASDIR}/PGC_ADHD_BP_OCD_MDD_SCZ3_PTSD/pgc-bip2021-all_noheader.vcf.tsv.gz \
--out ${DATADIR}/munged/Bipolar-Mullins_2021.parquet; fi


LABEL=Depression-Howard_2019
SUMSTATS=${GWASDIR}/PGC_ADHD_BP_OCD_MDD_SCZ3_PTSD/PGC_UKB_depression_genome-wide.txt.gz
SUMSTATS2=${GWASDIR}/PGC_ADHD_BP_OCD_MDD_SCZ3_PTSD/PGC_UKB_depression_genome-wide_rsid_chr.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then
zcat $SUMSTATS | cut -d ' ' -f1 > ${LABEL}_SNPs.txt 
zcat $SUMSTATS | tr ' ' '\t' > ${LABEL}.tmp.txt # extract SNPs to search
plink2 --memory 100000 --threads 4 --vcf ${VCFDIR}/All_20180423_b151_GRCh37p13.vcf.gz --extract ${LABEL}_SNPs.txt --make-just-pvar --out ${LABEL}
grep -Ev $'##' ${LABEL}.pvar | sed 's/#//g' | sed 's/ID/SNP/g' | grep -Ev $'X' | cut -f 1-3 > ${LABEL}.pvar2
cat <(head -n 1 ${LABEL}.pvar2 && tail -n +2 ${LABEL}.pvar2 | sort -k3,3 )> ${LABEL}.pvar3
cat <(head -n 1 ${LABEL}.tmp.txt && tail -n +2 ${LABEL}.tmp.txt | sort -k1,1 )> ${LABEL}.tmp2.txt
join -1 3 -2 1 --header ${LABEL}.pvar3 ${LABEL}.tmp2.txt | tr ' ' '\t' | gzip > ${SUMSTATS2}
python ${POLYFUNDIR}/munge_polyfun_sumstats.py --n 807553 --min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS2 --out ${DATADIR}/munged/${LABEL}.parquet; 
rm ${LABEL}.pvar ${LABEL}.pvar2 ${LABEL}.pvar3 ${LABEL}.tmp.txt ${LABEL}.tmp2.txt
fi


if [[ ! -f ${DATADIR}/munged/MDD-Wray_2018.parquet ]]; then
python ${POLYFUNDIR}/munge_polyfun_sumstats.py --min-info 0.6 --min-maf 0.001 \
--sumstats ${GWASDIR}/PGC_ADHD_BP_OCD_MDD_SCZ3_PTSD/PGC_MDD2018_ex23andMe_19fields.txt.gz \
--out ${DATADIR}/munged/MDD-Wray_2018.parquet
fi


if [[ ! -f ${DATADIR}/munged/OCD-IOCDF_2017.parquet ]]; then
python ${POLYFUNDIR}/munge_polyfun_sumstats.py --min-info 0.6 --min-maf 0.001 \
--n 9725 --sumstats ${GWASDIR}/PGC_ADHD_BP_OCD_MDD_SCZ3_PTSD/PGC_OCD_Aug2017/ocd_aug2017 \
--out ${DATADIR}/munged/OCD-IOCDF_2017.parquet
fi


if [[ ! -f ${DATADIR}/munged/PTSD_EUR-Nievergelt_2018.parquet ]]; then
zcat ${GWASDIR}/PGC_ADHD_BP_OCD_MDD_SCZ3_PTSD/pts_eur_freeze2_overall.results.gz | grep -Ev $'NA' | \
awk '!seen[$2]++' | gzip > ${GWASDIR}/PGC_ADHD_BP_OCD_MDD_SCZ3_PTSD/pts_eur_freeze2_overall.results.txt.gz 
python ${POLYFUNDIR}/munge_polyfun_sumstats.py --min-info 0.6 --min-maf 0.001 \
--sumstats ${GWASDIR}/PGC_ADHD_BP_OCD_MDD_SCZ3_PTSD/pts_eur_freeze2_overall.results.txt.gz \
--out ${DATADIR}/munged/PTSD_EUR-Nievergelt_2018.parquet
fi

LABEL=Schizophrenia-PGC3_2020
SUMSTATS=${GWASDIR}/PGC_ADHD_BP_OCD_MDD_SCZ3_PTSD/PGC3_SCZ_wave3_public.v2.tsv.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


#######################
## 3) Neuro traits NOS
LABEL=EduAttain-Lee_2018
SUMSTATS=${GWASDIR}/Intelligence_2018/GWAS_EA_excl23andMe.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--n 1131881 --min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


LABEL=Income-Hill_2019
SUMSTATS=${GWASDIR}/Household_Income_UKBiobank.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--n 286301 --min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


LABEL=BMI-Pulit_2019
SUMSTATS=${GWASDIR}/GIANT_traits/Bmi.giant-ukbb.meta-analysis.combined.23May2018.SNP_trimmed2.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then 
zcat ${GWASDIR}/GIANT_traits/Bmi.giant-ukbb.meta-analysis.combined.23May2018.SNP_trimmed.txt.gz | sed '/NA/d' |gzip > ${SUMSTATS}
python ${POLYFUNDIR}/munge_polyfun_sumstats.py --min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


LABEL=BrainVol-Jansen_2020
SUMSTATS=${GWASDIR}/brainvol/meta_analysis_BV_Jansenetal_2020.sumstats_noDup.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then 
zcat ${GWASDIR}/brainvol/meta_analysis_BV_Jansenetal_2020.sumstats.txt.gz | cut -f2 --complement |gzip > ${SUMSTATS}
python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
	--min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


LABEL=Neuroticism-Nagel_2018
SUMSTATS=${GWASDIR}/sumstats_neuroticism_ctg_format_noDup.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then 
zcat ${GWASDIR}/sumstats_neuroticism_ctg_format.txt.gz | cut -f2 --complement | sed '/NA/d' | gzip > ${SUMSTATS}
python ${POLYFUNDIR}/munge_polyfun_sumstats.py --min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; 
fi


LABEL=Intelligence-SavageJansen_2018
SUMSTATS=${GWASDIR}/Intelligence_2018/SavageJansen_2018_intelligence_metaanalysis.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
	--min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


LABEL=Worry-Nagel_2018
SUMSTATS=${GWASDIR}/sumstats_worry_ctg_format_noDup.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then
zcat ${GWASDIR}/sumstats_worry_ctg_format.txt.gz | cut -f2 --complement | gzip > ${SUMSTATS}
python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
	--min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


LABEL=Headache-Meng_2018
SUMSTATS=${GWASDIR}/headache2_2017-10-12_noN.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then
zcat ${GWASDIR}/headache2_2017-10-12.txt.gz | cut -d ' ' -f13 --complement | gzip > ${SUMSTATS}
python ${POLYFUNDIR}/munge_polyfun_sumstats.py --min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi



####################
## 4) Sleep traits 
LABEL=DaytimeSleepiness-Wang_2019
SUMSTATS=${GWASDIR}/Sleep_Genetics_GWAS/Saxena.fullUKBB.DaytimeSleepiness_adjBMI.sumstats.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--n 452071 --min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


LABEL=LongSleepDur-Dashti_2019
SUMSTATS=${GWASDIR}/Sleep_Genetics_GWAS/longsumstats.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--n 446118 --min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi

LABEL=ShortSleepDur-Dashti_2019
SUMSTATS=${GWASDIR}/Sleep_Genetics_GWAS/shortsumstats.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--n 446118 --min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


LABEL=Dozing-Jansen_2019
SUMSTATS=${GWASDIR}/Insomnia_Jansen_2019/Dozing_sumstats_Jansenetal2.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then 
zcat ${GWASDIR}/Insomnia_Jansen_2019/Dozing_sumstats_Jansenetal.txt.gz | sed '/X/d' | gzip > $SUMSTATS
python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


LABEL=Gettingup-Jansen_2019
SUMSTATS=${GWASDIR}/Insomnia_Jansen_2019/Gettingup_sumstats_Jansenetal2.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then 
zcat ${GWASDIR}/Insomnia_Jansen_2019/Gettingup_sumstats_Jansenetal.txt.gz | sed '/X/d' | gzip > $SUMSTATS
python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


LABEL=Insomnia-Jansen_2019
SUMSTATS=${GWASDIR}/Insomnia_Jansen_2019/Insomnia_sumstats_Jansenetal2.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then 
zcat ${GWASDIR}/Insomnia_Jansen_2019/Insomnia_sumstats_Jansenetal.txt.gz | sed '/X/d' | gzip > $SUMSTATS
python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


LABEL=Morningness-Jansen_2019
SUMSTATS=${GWASDIR}/Insomnia_Jansen_2019/Morningness_sumstats_Jansenetal2.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then
zcat ${GWASDIR}/Insomnia_Jansen_2019/Morningness_sumstats_Jansenetal.txt.gz | sed '/X/d' | gzip > $SUMSTATS
python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


LABEL=Napping-Jansen_2019
SUMSTATS=${GWASDIR}/Insomnia_Jansen_2019/Napping_sumstats_Jansenetal2.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then 
zcat ${GWASDIR}/Insomnia_Jansen_2019/Napping_sumstats_Jansenetal.txt.gz | sed '/X/d' | gzip > $SUMSTATS
python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


LABEL=Sleepdur-Jansen_2019
SUMSTATS=${GWASDIR}/Insomnia_Jansen_2019/Sleepdur_sumstats_Jansenetal2.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then 
zcat ${GWASDIR}/Insomnia_Jansen_2019/Sleepdur_sumstats_Jansenetal.txt.gz | sed '/X/d' | gzip > $SUMSTATS
python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


LABEL=Snoring-Jansen_2019
SUMSTATS=${GWASDIR}/Insomnia_Jansen_2019/Snoring_sumstats_Jansenetal2.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then 
zcat ${GWASDIR}/Insomnia_Jansen_2019/Snoring_sumstats_Jansenetal.txt.gz | sed '/X/d' | gzip > $SUMSTATS
python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi


###################
## 2) Psych traits 
LABEL=AlcoholDep-Walters_2018
SUMSTATS=${GWASDIR}/PGC_SUD_GWAS/PGC2020_AUDIT/pgc_alcdep.eur_discovery.aug2018_release.txt.gz
if [[ ! -f ${DATADIR}/munged/${LABEL}.parquet ]]; then 
python ${POLYFUNDIR}/munge_polyfun_sumstats.py \
--n 46568 --min-info 0.6 --min-maf 0.001 --sumstats $SUMSTATS --out ${DATADIR}/munged/${LABEL}.parquet; fi
