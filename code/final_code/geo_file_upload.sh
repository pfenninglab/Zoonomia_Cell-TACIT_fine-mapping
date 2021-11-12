FASTQDIR=/data/pfenninggroup/singleCell/Macaque_snATAC-seq/m015_Peanut/snATACseq
cd $FASTQDIR

## upload fastq folders
lftp ftp://geoftp:rebUzyi1@ftp-private.ncbi.nlm.nih.gov
cd uploads/badoi.phan@pitt.edu_01AX6gqO
mirror -R STA682A131

cd uploads/badoi.phan@pitt.edu_01AX6gqO
mirror -R STA682A134
exit


## upload ArchR project folders
PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
DATADIR=$PROJDIR/data/raw_data/rheMac10/Stauffer_caudate
cd $DATADIR

tar -czvf Stauffer_caudate_labeled_ArchR_Project.tar.gz ArchR_Stauffer_caudate_labeled

lftp ftp://geoftp:rebUzyi1@ftp-private.ncbi.nlm.nih.gov
cd uploads/badoi.phan@pitt.edu_01AX6gqO
mirror -R peak
put Stauffer_caudate_labeled_ArchR_Project.tar.gz



cd $PROJDIR/data/tidy_data/data_for_geo_upload
lftp ftp://geoftp:rebUzyi1@ftp-private.ncbi.nlm.nih.gov
cd uploads/badoi.phan@pitt.edu_01AX6gqO
put Stauffer_caudate_geo_submission_sheet.xlsx


