###############################
DIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/figures/explanatory/figure1_snATAC_cross-species_plots
mkdir -p $DIR/plots
cd $DIR/plots

# SNP blocks to plot
BLOCK1=$DIR/track_plots_mm10_regions.bed
TRACKS20=$DIR/step2_config_mouse_cSNAIL.ini

# tall track for figure
pyGenomeTracks --tracks $TRACKS20 --BED $BLOCK1 --fontSize 1 --width 2 --dpi 1000 --trackLabelFraction 0.1 --outFileName fig1_mouse_cSNAIL_pyGenomeTracks.pdf
pyGenomeTracks --tracks $TRACKS20 --BED $BLOCK1 --fontSize 1 --width 2 --dpi 1000 --trackLabelFraction 0.1 --outFileName fig1_mouse_cSNAIL_pyGenomeTracks40.pdf


