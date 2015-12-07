### Q1.glmnet.R created by Honglei Xie ##########################################################
### This script is to develop classification algorithm to predict which AML patients will have ####
### Complete Remission(CR) or will be Primary resistant(RESISTANT) using clinical data only #######

### NOTES #########################################################################################



### PREAMBLE ######################################################################################
# load required packages and scripts
library(randomForest);
library(BoutrosLab.statistics.classification);
library(BoutrosLab.prognosticsignature.general);
library(glmnet);
library(pROC);

source("~/cluster/svn/AlgorithmEvaluations/DREAM9/AML/GroupShare/maximize_bac.R");

### set woking dir and load training datasets #####################################################
setwd("~/cluster/svn/AlgorithmEvaluations/DREAM9/AML/Honglei/");

dir.data <- "~/cluster/svn/AlgorithmEvaluations/DREAM9/AML/GroupShare/Data/";
dir.data.test <- "~/isilon/private/Datasets/DREAM/DREAM9/AML/data/test_data";


train.data <- read.csv(file.path(dir.data, 'AML_processed.csv'));
scoringData <- read.csv(file.path(dir.data.test, 'scoringData-release.csv'));

### pre-processing datasets #######################################################################
# remove all outcome columns except for resp.simple which is the event of interest in sub-challenge 1
train.data <- train.data[, !names(train.data) %in% 
                           c('Overall_Survival', 'Relapse', 'Remission_Duration', 'vital.status')];
train.data <- train.data[, !names(train.data) %in% 
                             c('Relapse', 'Remission_Duration')];

# remove more columns
train.data <- train.data[, !names(train.data) %in% c('X.Patient_id')];
scoringData <- scoringData[, !names(scoringData) %in% c('X.Patient_id')];



# # checking missingness
# missing.prop <- colMeans(is.na(train.data));
# # using means to replace missing value
# mean.values <- colMeans(train.data[, which(missing.prop != 0)], na.rm = TRUE);
# 
# for (var.name in names(mean.values)) 
#   train.data[is.na(train.data[[var.name]]), var.name] <- mean.values[var.name];

# replace NAs with medians in scoringData
missing.prop <- colMeans(is.na(scoringData));
median.values <- apply(scoringData[, which(missing.prop != 0)], 
                       MARGIN = 2, 
                       FUN = function(x){median(x, na.rm = TRUE)});

for (var.name in names(median.values)) 
  scoringData[is.na(scoringData[[var.name]]), var.name] <- median.values[var.name];

# take subset of intersection of covariates in training dataset and scoring dataset
train.data <- train.data[
  train.data$cyto.cat %in% intersect(levels(train.data$cyto.cat), levels(scoringData$cyto.cat)) *
  train.data$Chemo.Simplest %in% intersect(levels(train.data$Chemo.Simplest) , levels(scoringData$Chemo.Simplest)) 
  == 1
  ,];

scoringData <- scoringData[
  scoringData$cyto.cat %in% intersect(levels(train.data$cyto.cat), levels(scoringData$cyto.cat)) *
  scoringData$Chemo.Simplest %in% intersect(levels(train.data$Chemo.Simplest) , levels(scoringData$Chemo.Simplest)) 
  == 1
  ,];

# process covariates 
train.data$cyto.cat <- factor(train.data$cyto.cat, 
                              levels = c('21','-5','-5,-7','-7','8','diploid','IM','inv16','inv9','Misc','t6;9','t8;21'));
train.data$Chemo.Simplest <- factor(train.data$Chemo.Simplest, 
                                    levels = c("Anthra-HDAC","Flu-HDAC","HDAC-Plus non Anthra","StdAraC-Plus"));

scoringData$cyto.cat <- factor(scoringData$cyto.cat, 
                               levels = c('21','-5','-5,-7','-7','8','diploid','IM','inv16','inv9','Misc','t6;9','t8;21'));

scoringData$Chemo.Simplest <- factor(scoringData$Chemo.Simplest, 
                                     levels = c("Anthra-HDAC","Flu-HDAC","HDAC-Plus non Anthra","StdAraC-Plus"));

# double check missingness, expecting both of them are all TRUE
table(colMeans(is.na(train.data)) == 0);
table(colMeans(is.na(scoringData)) == 0);


# only use clinical covariates
train.data.clin <- train.data[,c(1:36)];
scoringData.clin <- scoringData[,c(1:35)];
rownames(train.data.clin) <- NULL;
rownames(scoringData.clin) <- NULL;



### MODEL: fit a glm ##############################################################################
# randomly select 60 % of patients as training sets

performance <- data.frame(AUC = NA, max.BAC = NA);

for(i in 1:1000){
  i.train <- sample(1:nrow(train.data.clin), floor(nrow(train.data.clin)*0.6), replace = TRUE);
  nm.cov <- setdiff(names(train.data.clin), c('resp.simple'));
  
  design.matrix <- model.matrix( formula(paste('~', paste(nm.cov, collapse = '+'))), data = train.data.clin);
  
  obs.weights <- (1/prop.table(table(train.data.clin[i.train, 'resp.simple'])))[train.data.clin[i.train, 'resp.simple']];
  
  # fit a glm
  fit.glmnet <- cv.glmnet(
    x = design.matrix[i.train,], 
    y = as.numeric(train.data.clin[i.train, 'resp.simple'] == 'CR'), 
    weights = obs.weights, 
    family = 'binomial', 
    alpha = 1
  );
  
  # predict probability for the test dataset
  pred.glmnet.prob <- predict(fit.glmnet, newx = design.matrix[-i.train,], s = 'lambda.min', type = 'response');
  
  #true class for the training patients
  true.class <- train.data.clin$resp.simple[ i.train];
  true.class <- as.numeric(true.class == 'CR');
  
  model.roc <- roc(response = true.class, predictor = pred.glmnet.prob, auc = TRUE);
  performance[i,] <- rocplot(model.roc, filename = NULL)[1:2];
  
  
}
# calculate empirical 95% confidence interval
hist(performance$AUC, xlab = "AUC", main = " ");
hist(performance$max.BAC, xlab = "maxBAC", main = " ");

quantile(performance$AUC, probs = c(0.025, 0.975));
quantile(performance$max.BAC, probs = c(0.025, 0.975));

### apply for the scoring dataset #################################################################
design.matrix.sc <- model.matrix( formula(paste('~ ', paste(nm.cov, collapse = ' + '))),
                                  data = scoringData.clin);
prediction.prob <- predict(fit.glmnet, newx = design.matrix.sc, s = 'lambda.min', type = 'response');

# table the final predictions
prediction.class <- as.factor(apply(matrix(rep(NA, nrow(scoringData.clin) )), 
                              MARGIN = 2, 
                              FUN = function(x){ifelse(prediction.prob > 0.5, 'CR', 'RESISTENT')}
));

table(prediction.class);

### write session info file #######################################################################
save.session.profile(generate.filename('sub1_glmnet', 'Session.Profile', 'txt'));

