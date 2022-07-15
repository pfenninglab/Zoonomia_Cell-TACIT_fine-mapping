#!/bin/bash

PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/
DATADIR=${PROJDIR}/data/raw_data/reporter_assay/

SNP_BED=${DATADIR}/peaks/snp.bed

PHYLOP_OUT=${DATADIR}/peaks/snp_phyloP.bed
PHYLO=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/phyloP/human-centered-200m-Feb2021/200m_scoresPhyloP_20210214.11.bigWig

PHASTC_OUT=${DATADIR}/peaks/snp_phastcons.bed
PHAST=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/phyloP/Primates_PhastCons_scores_hg38/43prim_PhastCons.11.bigWig

bwtool extract bed ${SNP_BED} ${PHYLO} ${PHYLOP_OUT}
bwtool extract bed ${SNP_BED} ${PHAST} ${PHASTC_OUT}
