setwd("C:/Users/HCHO/Documents/MATLAB/keutsch-PlantHub/eng/scripts/R")
library(R.matlab)
library(car)
library(dplyr)
library(ggpubr)
require(lmodel2)
source('YorkFit.R')

##INITIALIZATION############################################
setwd("D:/Plant/RAW/190605.1")
data <- readMat("co_R_StepsAveraged.mat")

# Specify x and y columns and their associated standard deviations (SD)

x          <- data$co.blank        # Chamber out blank for each step
sigma_x    <- data$co.blank.std    # Standard deviation of the chamber out blank
y          <- data$co.plant        # Chamber out with plant for each step
sigma_y    <- data$co.plant.std    # Standard deviation of the chamber out with plant

leaf_area  <- data$co.leaf.area  # Leaf area
total_flow <- data$co.total.flow # Total flow

############################################################

plot(x,y)

print(paste("Number of Points (N): ", length(x), sep=""))

print("The York Fit Coefficients")
York_results = YorkFit(x,y,sigma_x,sigma_y,printCoefs=1, makeLine=1)

York_slope    <- York_results[2,1]
York_slope_se <- York_results[2,2]
York_int      <- York_results[1,1]
York_int_se   <- York_results[1,2]


# Compensation Point

comp_pt     <- York_int/(1-York_slope)
comp_pt_err <- (comp_pt)*sqrt((York_int_se/York_int)^2+(York_slope_se/York_slope)^2)
print(paste("Compensation Point: ",comp_pt, sep=''))
print(paste("CP Error: ",comp_pt_err,sep=''))

dep_vel <- (total_flow/leaf_area)*((1-York_slope)/York_slope)*(1/60)
print(paste("Deposition Velocity: ",dep_vel, sep=''))


# Calculate residuals
# First we need to generate bypass values for each chamber point given a step's bypass mean and std dev

yhat         <- York_slope*x + York_int  # Predicted values from fit

e_Yorkslopex <- sqrt((York_slope_se/York_slope)^2+(sigma_x/x)^2)*(York_slope)*x
e_yhat       <- sqrt(e_Yorkslopex^2 + York_int_se^2)

York_res     <- y - yhat
York_res_err <- sqrt(sigma_y^2 + e_yhat^2)
York_res_var <- York_res_err^2

plot(x,York_res_var)
title('Uniform Variance of Residuals - Visual Check')

## TESTING FOR NORMALITY
# Perform Shapiro-Wilk Test for testing normality of residuals
print(shapiro.test(York_res))

# Perform Kolmogrov-Smirnov Test to test for normality of residuals 
# print(ks.test(York_res,'pnorm',alternative = "two.sided"))

# Plot QQ Plot - QQ Plots help to show how the sample distribution compares to a normal (Gaussian) distribution
print(ggqqplot(York_res,title = 'Normality of the Residuals'))

# Technically, this applies to OLS only
print(ncvTest(lm(y~x)))


## TESTING FOR HOMOSCEDACITY (Uniform variance in both x and y)
# Want null hypothesis to show that the variances are the same; Use Levene's Test

#yhat         <- York_slope*x_rpts + York_int  # Predicted values from fit
#York_res     <- y_pts - yhat

#print(leveneTest(as.numeric(York_res),as.factor(ygroups),center=mean))

# Technically, this applies to OLS only
#print(ncvTest(lm(y~x)))

