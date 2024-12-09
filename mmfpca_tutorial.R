###########################################################################################
## Description: A step-by-step implementation of estimation algorithm for multilevel 
##              multivariate FPCA and associated procedures including generating simulation 
##              data and model fitting described in 'Joint Modeling of Evoked and Induced 
##              Event-Related Spectral Perturbations.'
###########################################################################################
## Main functions implemented:
## 1. eigenf_construct: Function constructing the multilevel multi- and uni-variate 
##                      two-dimensional eigenfunctions used in the simulation.
## 2. share_gen: Function that generates a multilevel multivariate functional data set
##               with shared eigenscores across different variates
## 3. corr_gen: Function that generates a multilevel multivariate functional data set
##              with eigenscores correlated between different variates
## 4. mmfpca: Function that fits multilevel multivariate FPCA for high-dimensional 
##            functional data
## 5. mufpca: Function that fits multilevel univaraite FPCA for high-dimensional functional
##            data
###########################################################################################
## Required files:
##    1. eigenfunction_construction.R
##    2. simulation_data_generation.R
##    3. mmfpca.R
##    4. face.Cov.mfpca.R  
##       available at https://github.com/refunders/refund/blob/master/R/face.Cov.mfpca.R
##    5. mmfpca_plots.R
###########################################################################################

###########################################################################################
# Set the working direction
###########################################################################################
setwd("path_to_file")          # replace "path_to_file" with the working direction

###########################################################################################
# Load packages and required files
###########################################################################################
require(Matrix)
require(MASS)
require(splines)
require(pracma)
require(mgcv)
require(mvtnorm)
require(refund)
require(ggplot2)
source("eigenfunction_construction.R")          # Functions that construct multi- and uni-variate eigenfunctions
source("simulation_data_generation.R")          # Functions that generate multilevel multivariate functional data sets
source("mmfpca.R")                              # Functions that fits multilevel multi- and uni-variate FPCA for high-dimensional functional data
source("face.Cov.mfpca.R")                      # Supporting functions for multilevel univariate FPCA can be downloaded
                                                # from https://github.com/refunders/refund/blob/master/R/face.Cov.mfpca.R
source("mmfpca_plots.R")                        # Functions that generate figures for the simulation results



###########################################################################################
# 1. Construct the multilevel multi- and univariate eigenfunctions
###########################################################################################
true_eigen = eigenf_construct()

###########################################################################################
# 2. One run of simulation with shared multivariate eigenscores
###########################################################################################

# 2.1 Generate a 50-subject 50-trial high-noise bivariate two-dimensional functional outcome 
#     using shared eigenscores #############################################################
shared_dt = shared_gen(eigen_lvl1_var1 = true_eigen$multi$lvl1$var1,
                       eigen_lvl1_var2 = true_eigen$multi$lvl1$var2,
                       eigen_lvl2_var1 = true_eigen$multi$lvl2$var1,
                       eigen_lvl2_var2 = true_eigen$multi$lvl2$var2,
                       z_mu_var1 = rep(0, 2500),
                       z_mu_var2 = rep(0, 2500),
                       U = 2^(-0.5*c(1:5) + 2.5),
                       V = 2^(-0.5*c(1:5) + 2),
                       sigma2_e = 1,
                       N = 50,
                       R = 50)

# 2.2 Fit the multilevel multi- and uni-variate FPCA ######################################
## multilevel multivariate FPCA
shared_multi = mmfpca(z_var1 = shared_dt$z_var1,
                      z_var2 = shared_dt$z_var2,
                      id_array = shared_dt$id_array,
                      x_axis = seq(0.02, 1, length = 50),
                      y_axis = seq(0.02, 1, length = 50),
                      mufpca_pve = 0.999,
                      mmfpca_pve = 0.95)

## multilevel univaraite FPCA for variate 1 and 2 separately
shared_uni_var1 = mufpca(z_matrix = shared_dt$z_var1,
                         id_array = shared_dt$id_array,
                         x_axis = seq(0.02, 1, length = 50),
                         y_axis = seq(0.02, 1, length = 50),
                         pve = 0.95)

shared_uni_var2 = mufpca(z_matrix = shared_dt$z_var2,
                         id_array = shared_dt$id_array,
                         x_axis = seq(0.02, 1, length = 50),
                         y_axis = seq(0.02, 1, length = 50),
                         pve = 0.95)

# 2.3 Compare the estimated eigenfunctions from the multilevel multi- and uni-variate FPCA
#     with the true multivariate eigenfunctions ###########################################
## subject-level
levelplot(t(matrix(true_eigen$multi$lvl1$var1[,1], 50)))
levelplot(t(matrix(shared_result$multi$eigen_lvl1_var1_multi_est[,1], 50)))
levelplot(t(matrix(shared_result$uni$eigen_lvl1_var1_uni_est[,1]/sqrt(2), 50)))
levelplot(t(matrix(true_eigen$multi$lvl1$var2[,1], 50)))
levelplot(t(matrix(shared_result$multi$eigen_lvl1_var2_multi_est[,1], 50)))
levelplot(t(matrix(shared_result$uni$eigen_lvl1_var2_uni_est[,1]/sqrt(2), 50)))

## trial-level
levelplot(t(matrix(true_eigen$multi$lvl2$var1[,1], 50)))
levelplot(t(matrix(shared_result$multi$eigen_lvl2_var1_multi_est[,1], 50)))
levelplot(t(matrix(shared_result$uni$eigen_lvl2_var1_uni_est[,1]/sqrt(2), 50)))
levelplot(t(matrix(true_eigen$multi$lvl2$var2[,1], 50)))
levelplot(t(matrix(shared_result$multi$eigen_lvl2_var2_multi_est[,1], 50)))
levelplot(t(matrix(shared_result$uni$eigen_lvl2_var2_uni_est[,1]/sqrt(2), 50)))

# 2.4 Compare the data reconstruction from the multilevel multi- and uni-variate FPCA ###
## randomly select a single trial from a subject
subject = sample(c(1:50), size = 1)
trial = sample(c(1:50), size = 1)

## variate 1
levelplot(t(matrix(shared_dt$z_var1[(50*(subject-1)+trial),], 50)))
levelplot(t(matrix(shared_result$multi$z_multi_pred_var1[(50*(subject-1)+trial),], 50)))
levelplot(t(matrix(shared_result$uni$z_uni_pred_var1[(50*(subject-1)+trial),], 50)))

## variate 2
levelplot(t(matrix(shared_dt$z_var2[(50*(subject-1)+trial),], 50)))
levelplot(t(matrix(shared_result$multi$z_multi_pred_var2[(50*(subject-1)+trial),], 50)))
levelplot(t(matrix(shared_result$uni$z_uni_pred_var2[(50*(subject-1)+trial),], 50)))



###########################################################################################
# 3. One run of simulation with shared multivariate eigenscores
###########################################################################################

# 3.1 Generate a 50-subject 50-trial high-noise bivariate two-dimensional functional outcome 
#     using correlated eigenscores ########################################################
corr_dt = corr_gen(eigen_lvl1_var1 = true_eigen$uni$lvl1$var1,
                   eigen_lvl1_var2 = true_eigen$uni$lvl1$var2,
                   eigen_lvl2_var1 = true_eigen$uni$lvl2$var1,
                   eigen_lvl2_var2 = true_eigen$uni$lvl2$var2,
                   z_mu_var1 = rep(0, 2500),
                   z_mu_var2 = rep(0, 2500),
                   U = 2^(-0.5*c(1:5) + 2.5),
                   V = 2^(-0.5*c(1:5) + 2),
                   rho_lvl1 = 0.5,
                   rho_lvl2 = 0.5,
                   sigma2_e = 1,
                   N = 50,
                   R = 50)

# 3.2 Fit the multilevel multi- and uni-variate FPCA ######################################
corr_result = mmfpca(z_var1 = corr_dt$z_var1,
                     z_var2 = corr_dt$z_var2,
                     id_array = corr_dt$id_array,
                     x_axis = seq(0.02, 1, length = 50),
                     y_axis = seq(0.02, 1, length = 50),
                     mufpca_pve = 0.999,
                     mmfpca_pve = 0.95)

# 2.3 Compare the estimated eigenfunctions from the multilevel multi- and uni-variate FPCA
#     with the true multivariate eigenfunctions ###########################################
## subject-level
levelplot(t(matrix(true_eigen$uni$lvl1$var1[,1], 50)))
levelplot(t(matrix(corr_result$multi$eigen_lvl1_var1_multi_est[,1]*sqrt(2), 50)))
levelplot(t(matrix(corr_result$uni$eigen_lvl1_var1_uni_est[,1], 50)))
levelplot(t(matrix(true_eigen$uni$lvl1$var2[,1], 50)))
levelplot(t(matrix(corr_result$multi$eigen_lvl1_var2_multi_est[,1]*sqrt(2), 50)))
levelplot(t(matrix(corr_result$uni$eigen_lvl1_var2_uni_est[,1], 50)))

## trial-level
levelplot(t(matrix(true_eigen$uni$lvl2$var1[,1], 50)))
levelplot(t(matrix(corr_result$multi$eigen_lvl2_var1_multi_est[,1]*sqrt(2), 50)))
levelplot(t(matrix(corr_result$uni$eigen_lvl2_var1_uni_est[,1], 50)))
levelplot(t(matrix(true_eigen$uni$lvl2$var2[,1], 50)))
levelplot(t(matrix(corr_result$multi$eigen_lvl2_var2_multi_est[,1]*sqrt(2), 50)))
levelplot(t(matrix(corr_result$uni$eigen_lvl2_var2_uni_est[,1], 50)))

# 2.4 Compare the data reconstruction from the multilevel multi- and uni-variate FPCA ###
## randomly select a single trial from a subject
subject = sample(c(1:50), size = 1)
trial = sample(c(1:50), size = 1)

## variate 1
levelplot(t(matrix(corr_dt$z_var1[(50*(subject-1)+trial),], 50)))
levelplot(t(matrix(corr_result$multi$z_multi_pred_var1[(50*(subject-1)+trial),], 50)))
levelplot(t(matrix(corr_result$uni$z_uni_pred_var1[(50*(subject-1)+trial),], 50)))

## variate 2
levelplot(t(matrix(corr_dt$z_var2[(50*(subject-1)+trial),], 50)))
levelplot(t(matrix(corr_result$multi$z_multi_pred_var2[(50*(subject-1)+trial),], 50)))
levelplot(t(matrix(corr_result$uni$z_uni_pred_var2[(50*(subject-1)+trial),], 50)))
