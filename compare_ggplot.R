#### plot ##########################################################################################
#### This script is to generate plots comparing model performance between three methods. ###########
####################################################################################################
library(ggplot2);
dataset <- data.frame(
    Group = c(rep("Downsampling BAC", train.times), rep("Downsampling AUC", train.times), 
              rep("Oversampling BAC", train.times), rep("Oversampling AUC", train.times)),
    Subplot = rep(c("BAC","AUC"), each = 1000, times = 2),
    Value = c(down_performance$BAC, down_performance$AUC, over_performance$BAC, over_performance$AUC)
    );

ggplot(dataset, 
    aes(x = Group, y = Value)) + geom_boxplot() +
    theme(text = element_text(size=20)) + 
    facet_wrap(~ Subplot, scales = "free_x")
