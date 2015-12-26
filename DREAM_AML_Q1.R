# SCRIPT NAME
# ===========
#   chipmunks_aml_q1.R
#
# PACKAGE REQUIRED
# ================
#   randomForest, pROC ,getopt
#
# MODEL SUMMARY
# =============
#   used all clinical features
#   re-categorized covariate `cyto.cat` with reference to published paper
#   used balanced random forest algorithm with parameterization tuned by cross validation
#   used both 10-fold CV and randomly splitting to validate results
#
# USAGE
# =====
#   For Unix-like system (tested) or Windows (not tested) run the following commands on the terminal:
#   cd DIRECTORY_CONTAIN_THIS_SCRIPT
#   Rscript --no-save --no-restore chipmunks_aml_q1.R -t PATH_TO_TRAINING_SET -s PATH_TO_TEST_SET -o DIRECTORY_TO_SAVE_OUTPUT
#   If in any case you don't have the `Rscript` command, you can run the following alternative command
#   R CMD BATCH --no-save --no-restore '--args -t PATH_TO_TRAINING_SET -s PATH_TO_TEST_SET -o DIRECTORY_TO_SAVE_OUTPUT' chipmunks_aml_q1.R
#
# INPUT/command line argument
# ===========================
#   --train  / -t (required) : path to the `csv` training data set (filename!)
#   --test   / -s (required) : path to the `csv` testing  data set (filename!)
#   --output / -o (optional) : an output `directory`(NOT filename!). Default to the current working directory 
#
# OUTPUT
# ======
#   two files (predictions on traning and testing data set are generated and saved in the provideded 
#       directory specified by output argument ('--output / -o'), with names:
#       Chipmunks-AML_subchallenge1_submission_train.csv
#       Chipmunks-AML_subchallenge1_submission_test.csv
###################################################################################################

### PREAMBLE ######################################################################################
# load required packages
library(randomForest);
library(pROC);
library(getopt);

### settings ######################################################################################
# get option from command line
params <- matrix(
    c(
        'train', 't', 1, 'character',
        'test' , 's', 1, 'character',
        'output', 'o', 2, 'character'
    	),
    ncol = 4,
    byrow = TRUE
    );

if (!interactive())  opt <- getopt(params);

# if NOT run from terminal / command line, these 3 file paths must be provided.
path.train.data <- opt$train;
path.test.data  <- opt$test;
path.output     <- if(is.null(opt$output)) '.' else opt$output;


### DEFINE AXILIARY FUNCTIONS #####################################################################

#### maxBAC.R ##################################################################
# DESCRIPTION
#   Find the threshold that maximizes the balanced accuracy
# INPUT
#   roc : a ROC object from pROC
# OUTPUT
#   a list of thresholds and max balanced accuracy value

maxBAC <- function(roc) {
    
    require(pROC);
    
    bac <- 0.5*roc$sensitivities + 0.5*roc$specificities;
    max.indices <- which.max(bac);
    out <- list(
        'threshold'  = roc$thresholds[max.indices], 
        'max.bac'    = max(bac),
        'coord'      = data.frame(
            'sensitivities' = roc$sensitivities[max.indices],
            'specificities' = roc$specificities[max.indices]
        	),
        'auc'        = pROC::auc(roc)
    	)
    class(out) = c('bac', class(out));
    return(out);
}


#### rocplot ##############################################################
# DESCRIPTION
#   Plot roc curve and compute AUC, balanced accuracy(BAC)
# INPUT
#   roc : a roc object from the pROC package
#   filename: a file(.tiff) to save the plot. If NULL, use the default graphic device
#   OTHERS : see pROC::plot.roc and graphics::plot for details. Recommend to use the dafaults
# OUTPUT
#   Generate a roc curves with point(s) giving best balanced accuracy highlighted. 
#   Also return a vector of AUC and max.BAC, the correspondent thresholds, sensitivities, specificities
#   Can use this to obtain maximum BAC and visualize the ROC at once

rocplot <- function(
    roc, filename = NULL, 
    print.thres = 'best', print.thres.col = 'red', 
    print.thres.cex = 2, print.thres.pattern.cex = 1.5, 
    grid = TRUE, cex.lab = 1.5, cex.axis = 1.5, 
    ...) {
    
    require(pROC);
    bac <- maxBAC(roc);
    
    if (!is.null(filename)) tiff(filename); 
    plot(
        roc,
        print.thres             = print.thres,
        print.thres.col         = print.thres.col,
        print.thres.cex         = print.thres.cex,
        print.thres.pattern.cex = print.thres.pattern.cex,
        grid                    = grid,
        cex.lab                 = cex.lab,
        cex.axis                = cex.axis,
        ...
    	);
    legend('bottomright', c(paste0('AUC: ', round(pROC::auc(roc), 3)), paste0('BAC: ', round(bac$max.bac, 3))), xjust = 1, cex = 1.5);
    if (!is.null(filename)) dev.off();
    
    return( data.frame(
        'AUC'         = pROC::auc(roc),
        'max.BAC'     = bac$max.bac,
        'threshold'   = bac$threshold,
        'sensitivity' = bac$coord$sensitivities,
        'specificity' = bac$coord$specificities
    		)	
    	);
}



### READ DATA and SOME PROCESSING  ################################################################
# remove all outcome columns except for "resp.simple" which is the event of interest in sub-challenge 1
train.data <- read.csv(path.train.data);

train.data <- train.data[, !names(train.data) %in% c('Overall_Survival', 'Relapse', 'Remission_Duration', 'vital.status')];
train.data <- train.data[, !names(train.data) %in% c('X.Patient_id')];

scoringData <- read.csv(path.test.data);

scoringData <- scoringData[, !names(scoringData) %in% c('X.Patient_id')];

# checking missingness
missing.prop <- colMeans(is.na(train.data));
# using means to replace missing value
mean.values <- colMeans(train.data[, which(missing.prop != 0)], na.rm = TRUE);

for (var.name in names(mean.values)) 
  train.data[is.na(train.data[[var.name]]), var.name] <- mean.values[var.name];

# replace NAs with medians in scoringData
missing.prop <- colMeans(is.na(scoringData));
median.values <- apply(scoringData[, which(missing.prop != 0)], 
                       MARGIN = 2, 
                       FUN = function(x){median(x, na.rm = TRUE)});

for (var.name in names(median.values)) 
    scoringData[is.na(scoringData[[var.name]]), var.name] <- median.values[var.name];

# double check missingness, expecting both of them are all TRUE
# table(colMeans(is.na(train.data)) == 0);
# table(colMeans(is.na(scoringData)) == 0);

# handle different levels of variables in training dataset and test dataset
train.data$Chemo.Simplest[which(train.data$Chemo.Simplest == "Anthra-Plus")] <- "Anthra-HDAC";
train.data$Chemo.Simplest <- factor(train.data$Chemo.Simplest, 
                                    levels = c("Anthra-HDAC", "Flu-HDAC", "HDAC-Plus non Anthra","StdAraC-Plus")
                                    );
# only use clinical features
train.data.clin <- train.data[,c(1:36)];
scoringData.clin <- scoringData[,c(1:35)];
rownames(train.data.clin) <- NULL;
rownames(scoringData.clin) <- NULL;

# visualization: run it to get the graph that has been inserted in the writeup
#spineplot(train.data$cyto.cat, train.data$resp.simple);

# recategorize "cyto.cat" based on the graph above
cate_cyto <- function(data){
    
    data$cyto.cat1 <- factor(data$cyto.cat %in% c('inv16','t8;21'));
    data$cyto.cat2 <- factor(data$cyto.cat %in% c('11q23','8','21','-5','-5,-7','-5,-7,+8','-7','-7,+8'));
    data$cyto.cat3 <- factor(data$cyto.cat %in% c('diploid','IM','inv9','Misc','t6;9','t9;22'));
    
    data$cyto.cat <- NULL;
    return(data)
	}

train.data.clin <- cate_cyto(train.data.clin);
scoringData.clin <- cate_cyto(scoringData.clin);



### MODEL: BALANCED RANDOM FOREST  ################################################################

### NOTE: ######################################################################################### 
# The validation and parameterization tuning process are omitted here.
###################################################################################################


# all clinical features 
var.selected <- colnames(scoringData.clin);
# response
resp.train <- as.factor(train.data[, 'resp.simple']);

# make sure you can produce the results exactly as same as what we have submitted
RANDOM_SEED <- 12345;
set.seed(RANDOM_SEED);

rf.fit <- randomForest(
    x = train.data.clin[, var.selected], 
    y = as.factor(train.data.clin[, 'resp.simple']),
    mtry = 2,
    ntree = 1e5,
    sampsize = rep(min(table(resp.train)), nlevels(resp.train))     
    );

# prediction for train dataset
submission_train <- predict(
    rf.fit, 
    newdata = train.data.clin, 
    type = 'prob'
	)[,1];

# update the best threshold (refer to writeup for more details)
train.true.class <- train.data.clin$resp.simple;
train.model.roc <- roc(response = train.true.class, predictor = submission_train, auc = TRUE);
threshold.new <- maxBAC(train.model.roc)$threshold;

# prediction for test dataset
rf.pred <- predict(
    rf.fit, 
    newdata = scoringData.clin, 
    type = 'prob'
	)[,1];

rf.pred  <- ifelse(rf.pred < threshold.new , 
                   (0.5/threshold.new)*rf.pred,
                   (0.5/(1- threshold.new))* rf.pred + 1 - (0.5/(1- threshold.new))
				);

### SAVE OUTPUTS ##################################################################################
outfile.train <- data.frame(
    '#Patient_id' = paste0("id_", sprintf("%03d", 1:nrow(train.data.clin))),
    'CR_Confidence' = submission_train,
    stringsAsFactors = FALSE,
    check.names = FALSE
	);

outfile.test <- data.frame(
    '#Patient_id' = paste0("id_", sprintf("%03d", 1:nrow(scoringData.clin))),
    'CR_Confidence' = rf.pred,
    stringsAsFactors = FALSE,
    check.names = FALSE
	);

write.csv(x = outfile.train, file = file.path(path.output, 'Chipmunks-AML_subchallenge1_submission_train.csv'), quote = FALSE, row.names = FALSE);
write.csv(x = outfile.test , file = file.path(path.output, 'Chipmunks-AML_subchallenge1_submission_test.csv' ), quote = FALSE, row.names = FALSE);

message(sprintf('Predctions on training dataset are saved as\n\t%s/Chipmunks-AML_subchallenge1_submission_train.csv', path.output));
message(sprintf('Predctions on testing dataset are saved as\n\t%s/Chipmunks-AML_subchallenge1_submission_test.csv', path.output));
message(paste0(rep('-', 100), collapse = ''));
message('DONE!!! CHEERS!!!\n');
