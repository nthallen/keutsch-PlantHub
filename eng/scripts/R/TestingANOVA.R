my_data <- PlantGrowth

my_data <- read.csv(file.choose())


library(dplyr)


# Compute the analysis of variance
res.aov <- aov(my_data$X0.3 ~ my_data$A, data = my_data)
# Summary of the analysis
summary(res.aov)
