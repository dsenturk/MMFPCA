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
##               with shared eigenscores across different variates (setting 1)
## 3. corr_gen: Function that generates a multilevel multivariate functional data set
##              with eigenscores correlated between different variates (setting 2)
## 4. mmfpca: Function that fits multilevel multivariate FPCA for high-dimensional 
##            functional data
## 5. mufpca: Function that fits multilevel univariate FPCA for high-dimensional functional
##            data
## 6. efunction_plots: Function that generates figures for the true eigenfunctions and the
##                     estimated eigenfunctions from the multilevel multi- and uni-variate
##                     FPCA
## 7. reconstruction_plots: Function that generates figures for the trial-specific functional
##                          data and reconstruction using multi-level multi- and uni-variate
##                          FPCA for a single trial
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
setwd("path_to_file")          # replace "path_to_file" with the path to working environment

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
require(tidyverse)
require(ggplot2)
source("eigenfunction_construction.R")          # Functions that construct multi- and uni-variate eigenfunctions
source("simulation_data_generation.R")          # Functions that generate multilevel multivariate functional data sets
source("mmfpca.R")                              # Functions that fits multilevel multi- and uni-variate FPCA for high-dimensional functional data
source("face.Cov.mfpca.R")                      # Supporting functions for multilevel univariate FPCA
                                                # available at https://github.com/refunders/refund/blob/master/R/face.Cov.mfpca.R
source("mmfpca_plots.R")                        # Functions that generate figures for the simulation results



###########################################################################################
# 1. Construct the multilevel multi- and univariate eigenfunctions
###########################################################################################
true_eigen = eigenf_construct()

###########################################################################################
# 2. One run of simulation under simulation setting 1
###########################################################################################

# 2.1 Generate a 50-subject 50-trial high-noise bivariate two-dimensional functional outcome 
#     under simulation setting 1 ###########################################################
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
                      mufpca_pve = 0.99,
                      mmfpca_pve = 0.95,
                      efunctions_multi = TRUE)

## multilevel univariate FPCA for variate 1 and 2 separately
shared_uni_var1 = mufpca(z_matrix = shared_dt$z_var1,
                         id_array = shared_dt$id_array,
                         x_axis = seq(0.02, 1, length = 50),
                         y_axis = seq(0.02, 1, length = 50),
                         pve = 0.95,
                         efunctions_multi = TRUE)

shared_uni_var2 = mufpca(z_matrix = shared_dt$z_var2,
                         id_array = shared_dt$id_array,
                         x_axis = seq(0.02, 1, length = 50),
                         y_axis = seq(0.02, 1, length = 50),
                         pve = 0.95,
                         efunctions_multi = TRUE)

# 2.3 Compare the estimated eigenfunctions from the multilevel multi- and uni-variate FPCA
#     with the true multivariate eigenfunctions ###########################################
## match the sign of estimated eigenfunctions to the true eigenfunctions
shared_efunctions_matched = efunctions_match(true_eigen$multi,
                                             shared_multi$efunctions,
                                             shared_uni_var1$efunctions,
                                             shared_uni_var2$efunctions)

## generate comparison plots of the true and estimated eigenfunctions
shared_figures = efunctions_plots(shared_efunctions_matched,
                                  x_axis = seq(0.02, 1, length = 50),
                                  y_axis = seq(0.02, 1, length = 50))

## subject level
print(shared_figures$lvl1$var1)              # Figure S3
print(shared_figures$lvl1$var2)              # Figure S4

## trial level
print(shared_figures$lvl2$var1)              # Figure S5
print(shared_figures$lvl2$var2)              # Figure S6


# 2.4 Compare the data reconstruction from the multilevel multi- and uni-variate FPCA ###
## randomly select a single trial from a subject
subject = sample(c(1:50), size = 1)
trial = sample(c(1:50), size = 1)

## generate the trial-specific data reconstruction figures for the selected trial
shared_reconstruction = reconstruction_plots(shared_dt,
                                             shared_multi,
                                             shared_uni_var1,
                                             shared_uni_var2,
                                             x_axis = seq(0.02, 1, length = 50),
                                             y_axis = seq(0.02, 1, length = 50),
                                             subject,
                                             trial)

print(shared_reconstruction)                 # Figure S7



###########################################################################################
# 3. One run of simulation under simulation setting 2
###########################################################################################

# 3.1 Generate a 50-subject 50-trial high-noise bivariate two-dimensional functional outcome 
#     under simulation setting 2 ###########################################################
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
## multilevel multivariate FPCA
corr_multi = mmfpca(z_var1 = corr_dt$z_var1,
                    z_var2 = corr_dt$z_var2,
                    id_array = corr_dt$id_array,
                    x_axis = seq(0.02, 1, length = 50),
                    y_axis = seq(0.02, 1, length = 50),
                    mufpca_pve = 0.99,
                    mmfpca_pve = 0.95,
                    efunctions_multi = FALSE)

## multilevel univariate FPCA for variate 1 and 2 separately
corr_uni_var1 = mufpca(z_matrix = corr_dt$z_var1,
                       id_array = corr_dt$id_array,
                       x_axis = seq(0.02, 1, length = 50),
                       y_axis = seq(0.02, 1, length = 50),
                       pve = 0.95,
                       efunctions_multi = FALSE)

corr_uni_var2 = mufpca(z_matrix = corr_dt$z_var2,
                       id_array = corr_dt$id_array,
                       x_axis = seq(0.02, 1, length = 50),
                       y_axis = seq(0.02, 1, length = 50),
                       pve = 0.95,
                       efunctions_multi = FALSE)

# 3.3 Compare the estimated eigenfunctions from the multilevel multi- and uni-variate FPCA
#     with the true multivariate eigenfunctions ###########################################
## match the sign of estimated eigenfunctions to the true eigenfunctions
corr_efunctions_matched = efunctions_match(true_eigen$uni,
                                           corr_multi$efunctions,
                                           corr_uni_var1$efunctions,
                                           corr_uni_var2$efunctions)

## generate comparison plots of the true and estimated eigenfunctions
corr_figures = efunctions_plots(corr_efunctions_matched,
                                x_axis = seq(0.02, 1, length = 50),
                                y_axis = seq(0.02, 1, length = 50))

## subject-level
print(corr_figures$lvl1$var1)              # Figure S8
print(corr_figures$lvl1$var2)              # Figure S9

## subject-level
print(corr_figures$lvl2$var1)              # Figure S10
print(corr_figures$lvl2$var2)              # Figure S11


# 3.4 Compare the data reconstruction from the multilevel multi- and uni-variate FPCA ###
## randomly select a single trial from a subject
subject = sample(c(1:50), size = 1)
trial = sample(c(1:50), size = 1)

## generate the trial-specific data reconstruction figures for the selected trial
corr_reconstruction = reconstruction_plots(corr_dt,
                                           corr_multi,
                                           corr_uni_var1,
                                           corr_uni_var2,
                                           x_axis = seq(0.02, 1, length = 50),
                                           y_axis = seq(0.02, 1, length = 50),
                                           subject,
                                           trial)

print(corr_reconstruction)                 # Figure S12
