Rscript --no-save --no-restore chipmunks_aml_q3.R -t ~/DREAM9/AML/data/training_data/trainingData-release.csv -s ~/DREAM9/AML/data/test_data/scoringData-release.csv 
# R CMD BATCH --no-save --no-restore '--args -t ~/DREAM9/AML/data/training_data/trainingData-release.csv -s ~/DREAM9/AML/data/test_data/scoringData-release.csv' chipmunks_aml_q3.R
