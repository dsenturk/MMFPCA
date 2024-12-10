###########################################################################################
## Description: Functions for generating figures for simulation results in 'Joint Modeling
##              of Evoked and Induced Event-Related Spectral Perturbations.'
###########################################################################################
## Functions included:
##    1. efunctions_match: Function that flips the estimated eigenfunctions to match with 
##                         the sign of the true eigenfunctions
##    2. efunctions_plot: Function that generates figures for the true eigenfunctions and the
##                        estimated eigenfunctions from the multilevel multi- and uni-variate
##                        FPCA
##    3. reconstruction_plots: Function that generates figures for a selected trial from a
##                             multilevel multivariate functional data set and its reconstruction
##                             using multi- and uni-variate FPCA
###########################################################################################



efunctions_match = function(
    efunctions_true,                    # the true subject- and trial-level bivariate eigenfunctions (list)
    efunctions_mmfpca,                  # the subject- and trial-level bivariate eigenfunctions estimated
                                        # from multilevel multivariate FPCA (list)
    efunctions_mufpca_var1,             # the subject- and trial-level univariate eigenfunctions estimated
                                        # from multilevel univariate FPCA for variate 1 (list)
    efunctions_mufpca_var2){            # the subject- and trial-level univariate eigenfunctions estimated
                                        # from multilevel univariate FPCA for variate 2 (list)

  #########################################################################################
  ## Description: Function that flips the estimated eigenfunctions to match with the sign
  ##              of the true eigenfunctions
  ## Args:        see above
  ## Returns:     true: a list containing the true eigenfunctions
  ##              multi: a list containing the eigenfunctions estimated from multilevel 
  ##                     multivariate FPCA
  ##              uni: a list containing the eigenfunctions estimated from multilevel
  ##                   univariate FPCA
  #########################################################################################    
  
  # denote the number of eigenfunctions at each level
  m1_true = ncol(efunctions_true$lvl1$var1)
  m2_true = ncol(efunctions_true$lvl2$var1)
  m1_multi = ncol(efunctions_mmfpca$lvl1$var1)
  m2_multi = ncol(efunctions_mmfpca$lvl2$var1)
  m1_uni_var1 = ncol(efunctions_mufpca_var1$lvl1)
  m2_uni_var1 = ncol(efunctions_mufpca_var1$lvl2)
  m1_uni_var2 = ncol(efunctions_mufpca_var2$lvl1)
  m2_uni_var2 = ncol(efunctions_mufpca_var2$lvl2)
  
  # match the sign of estimated eigenfunctions with the true eigenfunctions 
  # 1. Calculate mean squared deviation error (MSDE) of the estimated eigenfunctions 
  #    and the MSDE of the sign-flipped estimated eigenfunctions
  # 2. For each eigenfunction,
  #    If the MSDE of the original estimation is smaller, keep the original estimation
  #    Otherwise, flip the sign of the original estimation
  
  ## subject level
  for (m in 1:m1_true){
    # multivariate eigenfunctions
    if(m <= m1_multi){
      msde_og = sum((efunctions_true$lvl1$var1[,m] - efunctions_mmfpca$lvl1$var1[,m])^2) +
                sum((efunctions_true$lvl1$var2[,m] - efunctions_mmfpca$lvl1$var2[,m])^2) 
      msde_flip = sum((efunctions_true$lvl1$var1[,m] + efunctions_mmfpca$lvl1$var1[,m])^2) +
                  sum((efunctions_true$lvl1$var2[,m] + efunctions_mmfpca$lvl1$var2[,m])^2)
      if(msde_flip < msde_og){
        efunctions_mmfpca$lvl1$var1[,m] = - efunctions_mmfpca$lvl1$var1[,m]
        efunctions_mmfpca$lvl1$var2[,m] = - efunctions_mmfpca$lvl1$var2[,m]
      }
    }
    # univariate eigenfunctions for variate 1
    if(m <= m1_uni_var1){
      msde_og = sum((efunctions_true$lvl1$var1[,m] - efunctions_mufpca_var1$lvl1[,m])^2)
      msde_flip = sum((efunctions_true$lvl1$var1[,m] + efunctions_mufpca_var1$lvl1[,m])^2)
      if(msde_flip < msde_og){
        efunctions_mufpca_var1$lvl1[,m] = - efunctions_mufpca_var1$lvl1[,m]
      }
    }
    # univariate eigenfunctions for variate 2
    if(m <= m1_uni_var2){
      msde_og = sum((efunctions_true$lvl1$var2[,m] - efunctions_mufpca_var2$lvl1[,m])^2)
      msde_flip = sum((efunctions_true$lvl1$var2[,m] + efunctions_mufpca_var2$lvl1[,m])^2)
      if(msde_flip < msde_og){
        efunctions_mufpca_var2$lvl1[,m] = - efunctions_mufpca_var2$lvl1[,m]
      }
    }
  }
  
  ## trial level
  for (m in 1:m2_true){
    # multivariate eigenfunctions
    if(m <= m2_multi){
      msde_og = sum((efunctions_true$lvl2$var1[,m] - efunctions_mmfpca$lvl2$var1[,m])^2) +
        sum((efunctions_true$lvl2$var2[,m] - efunctions_mmfpca$lvl2$var2[,m])^2) 
      msde_flip = sum((efunctions_true$lvl2$var1[,m] + efunctions_mmfpca$lvl2$var1[,m])^2) +
        sum((efunctions_true$lvl2$var2[,m] + efunctions_mmfpca$lvl2$var2[,m])^2)
      if(msde_flip < msde_og){
        efunctions_mmfpca$lvl2$var1[,m] = - efunctions_mmfpca$lvl2$var1[,m]
        efunctions_mmfpca$lvl2$var2[,m] = - efunctions_mmfpca$lvl2$var2[,m]
      }
    }
    # univariate eigenfunctions for variate 1
    if(m <= m2_uni_var1){
      msde_og = sum((efunctions_true$lvl2$var1[,m] - efunctions_mufpca_var1$lvl2[,m])^2)
      msde_flip = sum((efunctions_true$lvl2$var1[,m] + efunctions_mufpca_var1$lvl2[,m])^2)
      if(msde_flip < msde_og){
        efunctions_mufpca_var1$lvl2[,m] = - efunctions_mufpca_var1$lvl2[,m]
      }
    }
    # univariate eigenfunctions for variate 2
    if(m <= m2_uni_var2){
      msde_og = sum((efunctions_true$lvl2$var2[,m] - efunctions_mufpca_var2$lvl2[,m])^2)
      msde_flip = sum((efunctions_true$lvl2$var2[,m] + efunctions_mufpca_var2$lvl2[,m])^2)
      if(msde_flip < msde_og){
        efunctions_mufpca_var2$lvl2[,m] = - efunctions_mufpca_var2$lvl2[,m]
      }
    }
  }
  
  # store results into a list
  uni = list()
  uni$lvl1$var1 = efunctions_mufpca_var1$lvl1
  uni$lvl1$var2 = efunctions_mufpca_var2$lvl1
  uni$lvl2$var1 = efunctions_mufpca_var1$lvl2
  uni$lvl2$var2 = efunctions_mufpca_var2$lvl2
  
  result = list(true = efunctions_true,
                multi = efunctions_mmfpca,
                uni = uni)
  
  return(result)
}



efunctions_plots = function(
    x_axis,                     # x axis of the two-dimensional functional domain (array of length M_x)
    y_axis,                     # y axis of the two-dimensional functional domain (array of length M_y)
    efunctions_matched){        # sign-matched true eigenfunctions and the estimation from multilevel
                                # multi- and uni-variate FPCA (list returned by 'efunctions_match')

  #########################################################################################
  ## Description: Function that generates figures for the true eigenfunctions and the
  ##              estimated eigenfunctions from the multilevel multi- and uni-variate
  ##              FPCA
  ## Definition:  M_x: number of sampling points on the x-axis of the two-dimensional functional domain
  ##              M_y: number of sampling points on the y-axis of the two-dimensional functional domain
  ## Args:        see above
  ## Returns:     a list of figures containing
  ##              true eigenfunctions and the estimated eigenfunctions from the multi- and
  ##              uni-variate FPCA at the subject and trial level for variate 1 and 2
  #########################################################################################      
  
  # denote the number of sampling points on each axis of the two-dimensional functional domain
  M_x = length(x_axis)
  M_y = length(y_axis)
  M = M_x*M_y
  
  # denote the number of eigenfunctions at each level
  m1_true = ncol(efunctions_matched$true$lvl1$var1)
  m2_true = ncol(efunctions_matched$true$lvl2$var1)
  m1_multi = ncol(efunctions_matched$multi$lvl1$var1)
  m2_multi = ncol(efunctions_matched$multi$lvl2$var1)
  m1_uni_var1 = ncol(efunctions_matched$un$lvl1$var1)
  m2_uni_var1 = ncol(efunctions_matched$un$lvl2$var1)
  m1_uni_var2 = ncol(efunctions_matched$un$lvl1$var2)
  m2_uni_var2 = ncol(efunctions_matched$un$lvl2$var2)
  
  # generate plots comparing the true eigenfunctions and the multi- and uni-variate estimation at each level for each variate
  
  ## prepare data for plotting the subject-level eigenfunctions and estimations
  ### subject level: variate 1
  efunctions_lvl1_var1 = data.frame(z = c(efunctions_matched$true$lvl1$var1),
                                    x = rep(rep(x_axis, each = M_y), m1_true),
                                    y = rep(rep(y_axis, M_x), m1_true),
                                    number = rep(c("1st", "2nd", "3nd", paste0(4:m1_true, "th")),
                                                   each = M),
                                    type = "True") %>%
    union(data.frame(z = c(efunctions_matched$multi$lvl1$var1[,1:min(m1_true, m1_multi)]),
                     x = rep(rep(x_axis, each = M_y), min(m1_true, m1_multi)),
                     y = rep(rep(y_axis, M_x), min(m1_true, m1_multi)),
                     number = rep(c("1st", "2nd", "3nd",paste0(4:min(m1_true, m1_multi), "th")),
                                  each = M),
                     type = "Multivariate estimation")) %>%
    union(data.frame(z = c(efunctions_matched$uni$lvl1$var1[,1:min(m1_true, m1_uni_var1)]),
                     x = rep(rep(x_axis, each = M_y), min(m1_true, m1_uni_var1)),
                     y = rep(rep(y_axis, M_x), min(m1_true, m1_uni_var1)),
                     number = rep(c("1st", "2nd", "3nd", paste0(4:min(m1_true, m1_uni_var1), "th")),
                                  each = M),
                     type = "Univariate estimation"))
  ### subject level: variate 2
  efunctions_lvl1_var2 = data.frame(z = c(efunctions_matched$true$lvl1$var2),
                                    x = rep(rep(x_axis, each = M_y), m1_true),
                                    y = rep(rep(y_axis, M_x), m1_true),
                                    number = rep(c("1st", "2nd", "3nd", paste0(4:m1_true, "th")),
                                                 each = M),
                                    type = "True") %>%
    union(data.frame(z = c(efunctions_matched$multi$lvl1$var2[,1:min(m1_true, m1_multi)]),
                     x = rep(rep(x_axis, each = M_y), min(m1_true, m1_multi)),
                     y = rep(rep(y_axis, M_x), min(m1_true, m1_multi)),
                     number = rep(c("1st", "2nd", "3nd",paste0(4:min(m1_true, m1_multi), "th")),
                                  each = M),
                     type = "Multivariate estimation")) %>%
    union(data.frame(z = c(efunctions_matched$uni$lvl1$var2[,1:min(m1_true, m1_uni_var2)]),
                     x = rep(rep(x_axis, each = M_y), min(m1_true, m1_uni_var2)),
                     y = rep(rep(y_axis, M_x), min(m1_true, m1_uni_var2)),
                     number = rep(c("1st", "2nd", "3nd", paste0(4:min(m1_true, m1_uni_var2), "th")),
                                  each = M),
                     type = "Univariate estimation"))
  ### trial level: variate 1
  efunctions_lvl2_var1 = data.frame(z = c(efunctions_matched$true$lvl2$var1),
                                    x = rep(rep(x_axis, each = M_y), m2_true),
                                    y = rep(rep(y_axis, M_x), m2_true),
                                    number = rep(c("1st", "2nd", "3nd", paste0(4:m2_true, "th")),
                                                 each = M),
                                    type = "True") %>%
    union(data.frame(z = c(efunctions_matched$multi$lvl2$var1[,1:min(m2_true, m2_multi)]),
                     x = rep(rep(x_axis, each = M_y), min(m2_true, m2_multi)),
                     y = rep(rep(y_axis, M_x), min(m2_true, m2_multi)),
                     number = rep(c("1st", "2nd", "3nd",paste0(4:min(m2_true, m2_multi), "th")),
                                  each = M),
                     type = "Multivariate estimation")) %>%
    union(data.frame(z = c(efunctions_matched$uni$lvl2$var1[,1:min(m2_true, m2_uni_var1)]),
                     x = rep(rep(x_axis, each = M_y), min(m2_true, m2_uni_var1)),
                     y = rep(rep(y_axis, M_x), min(m2_true, m2_uni_var1)),
                     number = rep(c("1st", "2nd", "3nd", paste0(4:min(m2_true, m2_uni_var1), "th")),
                                  each = M),
                     type = "Univariate estimation"))
  ### trial level: variate 2
  efunctions_lvl2_var2 = data.frame(z = c(efunctions_matched$true$lvl2$var2),
                                    x = rep(rep(x_axis, each = M_y), m2_true),
                                    y = rep(rep(y_axis, M_x), m2_true),
                                    number = rep(c("1st", "2nd", "3nd", paste0(4:m2_true, "th")),
                                                 each = M),
                                    type = "True") %>%
    union(data.frame(z = c(efunctions_matched$multi$lvl2$var2[,1:min(m2_true, m2_multi)]),
                     x = rep(rep(x_axis, each = M_y), min(m2_true, m2_multi)),
                     y = rep(rep(y_axis, M_x), min(m2_true, m2_multi)),
                     number = rep(c("1st", "2nd", "3nd",paste0(4:min(m2_true, m2_multi), "th")),
                                  each = M),
                     type = "Multivariate estimation")) %>%
    union(data.frame(z = c(efunctions_matched$uni$lvl2$var2[,1:min(m2_true, m2_uni_var2)]),
                     x = rep(rep(x_axis, each = M_y), min(m2_true, m2_uni_var2)),
                     y = rep(rep(y_axis, M_x), min(m2_true, m2_uni_var2)),
                     number = rep(c("1st", "2nd", "3nd", paste0(4:min(m2_true, m2_uni_var2), "th")),
                                  each = M),
                     type = "Univariate estimation"))
  
  ## ensure 'type' is a factor with the desired level order
  efunctions_lvl1_var1$type = factor(efunctions_lvl1_var1$type, 
                                      levels = c("True", "Multivariate estimation", "Univariate estimation"))
  efunctions_lvl1_var2$type = factor(efunctions_lvl1_var2$type, 
                                     levels = c("True", "Multivariate estimation", "Univariate estimation"))
  efunctions_lvl2_var1$type = factor(efunctions_lvl2_var1$type, 
                                     levels = c("True", "Multivariate estimation", "Univariate estimation"))
  efunctions_lvl2_var2$type = factor(efunctions_lvl2_var2$type, 
                                     levels = c("True", "Multivariate estimation", "Univariate estimation"))
  
  ## generate the plots 
  fig_lvl1_var1 = ggplot(efunctions_lvl1_var1) +
    geom_tile(aes(x = x, y = y, fill = z)) +
    facet_grid(rows = vars(number), cols = vars(type)) +
    scale_fill_distiller(palette = "RdBu") +
    labs(x = "", y = "", fill = "", title = "subject-level eigenfunctions of variate 1") +
    theme_minimal()
  fig_lvl1_var2 = ggplot(efunctions_lvl1_var2) +
    geom_tile(aes(x = x, y = y, fill = z)) +
    facet_grid(rows = vars(number), cols = vars(type)) +
    scale_fill_distiller(palette = "RdBu") +
    labs(x = "", y = "", fill = "", title = "subject-level eigenfunctions of variate 2") +
    theme_minimal()
  fig_lvl2_var1 = ggplot(efunctions_lvl1_var1) +
    geom_tile(aes(x = x, y = y, fill = z)) +
    facet_grid(rows = vars(number), cols = vars(type)) +
    scale_fill_distiller(palette = "RdBu") +
    labs(x = "", y = "", fill = "", title = "trial-level eigenfunctions of variate 1") +
    theme_minimal()
  fig_lvl2_var2 = ggplot(efunctions_lvl2_var2) +
    geom_tile(aes(x = x, y = y, fill = z)) +
    facet_grid(rows = vars(number), cols = vars(type)) +
    scale_fill_distiller(palette = "RdBu") +
    labs(x = "", y = "", fill = "", title = "trial-level eigenfunctions of variate 2") +
    theme_minimal()
  
  # store all figures into a list
  fig_list = list()
  fig_list$lvl1$var1 = fig_lvl1_var1
  fig_list$lvl1$var2 = fig_lvl1_var2
  fig_list$lvl2$var1 = fig_lvl2_var1
  fig_list$lvl2$var2 = fig_lvl2_var2
  
  return(fig_list)
  
}


reconstruction_plots = function(
    z_data,                     # a multilevel multivariate two-dimensional functional data set generated
                                # by 'shared_gen' or 'corr_gen' (list)
    mod_mmfpca,                 # result of fitting multilevel multivariate FPCA for the multilevel multivariate
                                # functional data (z_data) using 'mmfpca' (list)
    mod_mufpca_var1,            # result of fitting multilevel univariate FPCA for the first variate of the 
                                # multilevel multivariate functional data (z_data) using 'mufpca' (list)
    mod_mufpca_var2,            # result of fitting multilevel univariate FPCA for the second variate of the 
                                # multilevel multivariate functional data (z_data) using 'mufpca' (list)
    x_axis,                     # x axis of the two-dimensional functional domain (array of length M_x)
    y_axis,                     # y axis of the two-dimensional functional domain (array of length M_y)
    subject,                    # index for the selected subject (scalar)
    trial){                     # index for the selected trial within the selected subject (scalar)

  #########################################################################################
  ## Description: Function that generates figures for a selected trial from a multilevel 
  ##              multivariate functional data set and its reconstruction using multi- and 
  ##              uni-variate FPCA
  ## Definition:  M_x: number of sampling points on the x-axis of the two-dimensional functional domain
  ##              M_y: number of sampling points on the y-axis of the two-dimensional functional domain
  ## Args:        see above
  ## Returns:     a figure containing
  ##              a selected trial from a multilevel multivariate functional data set and 
  ##              its reconstruction using multi- and uni-variate FPCA
  #########################################################################################      
  
  # denote the number of sampling points on each axis of the two-dimensional functional domain
  M_x = length(x_axis)
  M_y = length(y_axis)
  M = M_x*M_y
  # locate the selected trial in the data
  N = length(unique(z_data$id_array))                  # number of subjects
  R_list = array(dim = N)                              # number of trials for each subject
  for (i in 1:N){
    R_list[i] = sum(as.numeric(id_array == i))
  }
  if(subject == 1){
    trial_index = trial
  }else{
    trial_index = trial + sum(R_list[1:(subject-1)])
  }
  
  # prepare data for plotting the trial-specific data and reconstructions
  trial_dt = data.frame(z = c(z_data$z_var1[trial_index,]),
                        x = rep(x_axis, each = M_y),
                        y = rep(y_axis, M_x),
                        variate = "Variate 1",
                        type = "Raw data") %>%
    union(data.frame(z = c(z_data$z_var2[trial_index,]),
                     x = rep(x_axis, each = M_y),
                     y = rep(y_axis, M_x),
                     variate = "Variate 2",
                     type = "Raw data")) %>%
    union(data.frame(z = c(mod_mmfpca$z_hat$var1[trial_index,]),
                     x = rep(x_axis, each = M_y),
                     y = rep(y_axis, M_x),
                     variate = "Variate 1",
                     type = "Multivariate reconstruction")) %>%
    union(data.frame(z = c(mod_mmfpca$z_hat$var2[trial_index,]),
                     x = rep(x_axis, each = M_y),
                     y = rep(y_axis, M_x),
                     variate = "Variate 2",
                     type = "Multivariate reconstruction")) %>%
    union(data.frame(z = c(mod_mufpca_var1$z_hat[trial_index,]),
                     x = rep(x_axis, each = M_y),
                     y = rep(y_axis, M_x),
                     variate = "Variate 1",
                     type = "Univariate reconstruction")) %>%
    union(data.frame(z = c(mod_mufpca_var2$z_hat[trial_index,]),
                     x = rep(x_axis, each = M_y),
                     y = rep(y_axis, M_x),
                     variate = "Variate 2",
                     type = "Univariate reconstruction"))
  trial_dt$type = factor(trial_dt$type,
                         levels = c("Raw data", "Multivariate reconstruction", "Univariate reconstruction"))
  
  # generate the plots
  fig_trial = ggplot(trial_dt) +
    geom_tile(aes(x = x, y = y, fill = z)) +
    facet_grid(rows = vars(type), cols = vars(variate)) +
    scale_fill_distiller(palette = "RdBu") +
    labs(x = "", y = "", fill = "", title = paste("trial", trial, "of subject", subject)) +
    theme_minimal()
  
  return(fig_trial)
  
}

