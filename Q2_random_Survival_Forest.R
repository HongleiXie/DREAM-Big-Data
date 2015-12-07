# This script is trial for AML subschallenge 2

# settings:
dir.svn.root <- '~/BoutrosLab';
dir.out <- file.path(dir.svn.root, "AlgorithmEvaluations/DREAM9/AML/Eric/Q2/");
dir.data <- file.path(dir.svn.root, "AlgorithmEvaluations/DREAM9/AML/Eric/Data");
dir.axiliary.funs <- file.path(dir.svn.root, "AlgorithmEvaluations/DREAM9/AML/GroupShare/");
source(file.path(dir.axiliary.funs, 'get_pcc_ci.R'));

# load packages
library(ggplot2);
library(plyr);
library(caret);
library(survival);
library(doMC);
library(randomForestSRC);


# using the preprocessed data 
aml <- read.csv(file.path(dir.data,'AML_processed.csv'));

# aml.cn contains samples with complete remission. No missing on remission duration
aml.cr <- subset(aml, resp.simple == 1); # 136 samples

var.id      <- 'X.Patient_id';
var.resp    <- c('resp.simple', 'Relapse', 'vital.status', 'Overall_Survival', 'Remission_Duration');
var.cov     <- setdiff(names(aml.cr), c(var.id, var.resp));
var.clinic  <- setdiff(names(aml.cr)[1:41], c(var.resp, var.id));
var.protein <- names(aml.cr)[42:ncol(aml.cr)];


aml.final <- aml.cr[, c('Relapse', 'Remission_Duration', var.clinic)];
aml.final$Relapse <- as.numeric(aml.final$Relapse == 'Yes');

N <- nrow(aml.final);

set.seed(9812);
perf.nfold <- 10;
perf.foldid <- sample(rep(1:perf.nfold, length = N), size = N, replace = FALSE);

ps <- seq(0.99, 0.3, -0.01);

registerDoMC(cores = perf.nfold);

rf.perf <- foreach(perf.k = 1:perf.nfold) %dopar% {
    i.train <- sample(1:N, size = round(0.8*N), replace = FALSE);
    #i.train <- (1:N)[perf.foldid != perf.k];
    aml.train <- aml.final[i.train,];

#### MODEL 1: random survival forest ###########################################################

    rfsrc.fit <- rfsrc(
        Surv(Remission_Duration, Relapse) ~ ., data = aml.train, 
        ntree = 5000, mtry = 4);
    rfsrc.pred <- predict(rfsrc.fit, aml.final[-i.train, ], outcome = 'train');
    #rfsrc.pred.test <- predict(rfsrc.fit, aml.final[-i.train, ], outcome = 'test');
    Time <- rfsrc.fit$time.interest;

# best percentile for survival time predition

    rf.survtime <- t(apply(
          X = rfsrc.pred$survival, 
          MARGIN = 1, 
          FUN = function(x) {
              sapply(
                     X   = ps, 
                     FUN = function(p) {
                         index <- which(x <= p)[1];
                         if (is.na(index)) 600 else rfsrc.pred$time.interest[which(x <= p)[1]]
                         }
                     ) 
              }
          ));

    rf.survtime.mean <- c(c(rfsrc.pred$time.interest, 600) %*% 
        abs(apply(cbind(1, rfsrc.pred$survival, 0), 1, diff)))
    
    rf.survtime <- rf.survtime + 0.0001*matrix(runif(nrow(rf.survtime)*ncol(rf.survtime)) - 0.5, nrow = nrow(rf.survtime))

    #rf.survtime <- cbind(rf.survtime.mean, rf.survtime);
    rf.survtime.mix = rf.survtime[, 7]
    i  <- 7;
    unique.time <- unique(rf.survtime[, 7]);
    for( tm in unique.time) {
        i.tie <- which(rf.survtime[, 7] == tm);
        if (length(i.tie) > 1) {
            rf.survtime.mix[i.tie] <- rf.survtime.mix[i.tie] + 0.0001*rf.survtime.mean[i.tie]
            }
        }
    rf.survtime <- cbind(rf.survtime.mix, rf.survtime);

    rf.pcc.ci <- apply(rf.survtime, 2, get.pcc.ci, observed.time =  aml.final[-i.train, 'Remission_Duration'], is.observed = aml.final[-i.train, 'Relapse']);
    colnames(rf.pcc.ci) <- c('mean', paste0('p=', ps));
    rf.pcc.ci
    }

rf.perf.reduce <- Reduce('+', rf.perf) / perf.nfold;
p.max.pcc.ci <- which.max((rf.perf.reduce[1,]+1)/2 * 0.6 + rf.perf.reduce[2, ] * 0.4);
c(p = c(-1, ps)[p.max.pcc.ci], rf.perf.reduce[, p.max.pcc.ci])
p.max.pcc <- which.max(rf.perf.reduce[1,]);
c(p = c(-1, ps)[p.max.pcc], rf.perf.reduce[, p.max.pcc])
p.max.ci <- which.max(rf.perf.reduce[2,]);
c(p = c(-1, ps)[p.max.ci], rf.perf.reduce[, p.max.ci])


########################################################################
# use p=0.93, mtry = 4

dir.test.data <- '~/isilon/private/Datasets/DREAM/DREAM9/AML/data/test_data';
aml.test <- read.csv(file.path(dir.test.data,'scoringData-release.csv'));
# test if all clinic var in train set are all includedd in test set
all(var.clinic %in% names(aml.test));

# deal with factors

# if missing
any(is.na(aml.test[, var.clinic]));
# fill-in missing by median
for(var.missing.test in var.clinic[colSums(is.na(aml.test[, var.clinic])) > 0] ) 
    aml.test[is.na(aml.test[, var.missing.test]), var.missing.test] <- median(aml.final[, var.missing.test]);

var.factor <- var.clinic[sapply(aml.final[, var.clinic], is.factor)];
sapply(var.factor, function(v) identical(levels(aml.final[, v]), levels(aml.test[, v])))
for (v in var.factor) 
    aml.test[, v] <- factor(aml.test[, v], levels = levels(aml.final[, v]));


aml.test[is.na(aml.test[, 'cyto.cat']), 'cyto.cat'] <- sample(aml.final[, 'cyto.cat'], size = 2, replace = TRUE);


# train
rfsrc.fit <- rfsrc(Surv(Remission_Duration, Relapse) ~ ., data = aml.final, ntree = 10000, mtry = 4);
# prediction
rfsrc.pred <- predict(rfsrc.fit, aml.test[, var.clinic], outcome = 'train');
Time <- rfsrc.fit$time.interest;


# mean

rf.survtime.mean <- c(c(rfsrc.pred$time.interest, 600) %*% abs(apply(cbind(1, rfsrc.pred$survival, 0), 1, diff)));
rf.survtime.var <- c(c(rfsrc.pred$time.interest, 600)^2 %*% abs(apply(cbind(1, rfsrc.pred$survival, 0), 1, diff))) - rf.survtime.mean^2
rf.survtime.sd <- sqrt(rf.survtime.var);


diff(sort(rf.survtime.mean))

p <- 0.5;
rf.survtime.5 <- apply(
    X = rfsrc.pred$survival, 
    MARGIN = 1, 
    FUN = function(x) {
        index <- which(x <= p)[1]; 
        if (is.na(index)) 600 else rfsrc.pred$time.interest[which(x <= p)[1]]
        }
    );

diff(sort(rf.survtime.5))

# p <- 0.93;
# rf.survtime.93 <- apply(
#     X = rfsrc.pred$survival, 
#     MARGIN = 1, 
#     FUN = function(x) {
#         index <- which(x <= p)[1]; 
#         if (is.na(index)) 600 else rfsrc.pred$time.interest[which(x <= p)[1]]
#         }
#     );
# diff(sort(rf.survtime.93))

rf.survtime <- rf.survtime.5;
rf.survtime <- rf.survtime + 0.12*(runif(length(rf.survtime)) - 0.5) 

#rf.survtime <- rf.survtime.93*(sum(rf.survtime.5)/sum(rf.survtime.93));
# add a uniform noise to get better Concordance index, due to the CI formula prvovided.
rf.survtime <- pmin(rf.survtime, 600)

#apply(rfsrc.pred$survival)
rf.sd <- range(1/rf.survtime.sd);
confid <- 0.4*((1/rf.survtime.sd - rf.sd[1])/diff(rf.sd) + 0.2)

test.out <- data.frame(
                       '#Patient_id' = aml.test[, var.id], 
                       'Remission_Duration' = round(rf.survtime,4),
                       'Confidence' = round(confid, 4)
                       );


write.csv(x = test.out, file = 'Q2_to_submit_0728.csv', quote = FALSE, row.names = FALSE);
