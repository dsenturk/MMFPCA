###########################################################################################
## Description: Functions for fitting multilevel multivariate FPCA via estimation algorithm
##              described in 'Joint Modeling of Evoked and Induced Event-Related Spectral 
##              Perturbations.'
###########################################################################################
## Functions included:
## Main functions:
##    1. mmfpca: Function that fits multilevel multivariate FPCA for two-dimensional 
##               functional data
##    2. mufpca: Function that fits multilevel univariate FPCA using fast covariance  
##               estimation (FACE) method tailored for two-dimensional functional data 
## Supporting functions used by main function:
##    1. mfpca.face_center: Function that fits multilevel univariate FPCA for one-dimensional
##                          functional observations using FACE method (Cui et al. 2023) given
##                          pre-estimated mean function
##       * this function is modified from the 'mfpca.face' function in the 'refund' package
##       * the overall mean function needs to be pre-estimated and provided in the input
##         (different from the original 'mfpca.face' function)
##       * supporting functions utilized in the 'mfpca.face_center' function can be downloaded 
##         from https://github.com/refunders/refund/blob/master/R/face.Cov.mfpca.R and
##         manually loaded into the environment via 'source('face.Cov.mfpca.R')'
###########################################################################################



mmfpca = function(
    z_var1,                 # multilevel (vectorized two-dimensional) functional data from variate 1 
                            # (matrix of dimension M*R_sum)
    z_var2,                 # multilevel (vectorized two-dimensional) functional data from variate 2 
                            # (matrix of dimension M*R_sum)
    id_array,               # id information for each trial to identify the subject it belong to (array of length R_sum)
    x_axis,                 # x axis of the two-dimensional functional domain (array of length M_x)
    y_axis,                 # y axis of the two-dimensional functional domain (array of length M_y)
    mufpca_pve,             # percentage of variance explained (pve) threshold for multilevel univariate FPCA (scalar < 1)
    mmfpca_pve,             # pve threshold for multilevel multivariate FPCA (scalar < 1)
    efunctions_multi = T){  # indicator of whether the data is generated using multivariate eigenfunctions 
                            # (logical, TRUE = using multivariate eigenfunctions, FALSE = using univariate eigenfunctions)
  
  #########################################################################################
  ## Description: Function that fits multilevel multivariate FPCA for two-dimensional functional data
  ## Definition:  M_x: number of sampling points on the x-axis of the two-dimensional functional domain
  ##              M_y: number of sampling points on the y-axis of the two-dimensional functional domain
  ##              M = M_x * M_y: total number of sampling points on the two-dimensional functional domain
  ##              N: number of subjects
  ##              R_sum: total number of trials
  ## Args:        see above
  ## Returns:     mod_mmfpca: a list containing the estimated multilevel multivariate 
  ##                          eigencomponents and the corresponding data reconstruction
  #########################################################################################  
  
  N = length(unique(id_array))                         # number of subjects
  R_list = array(dim = N)                              # number of trials for each subject
  for (i in 1:N){
    R_list[i] = sum(as.numeric(id_array == i))
  }
  R_sum = length(id_array)                             # total number of trials
  M = dim(z_var1)[2]                                   # total number of sampling points
  
  # estimate the smoothed two-dimensional mean function
  ## calculate the empirical mean
  z_var1_mean = apply(z_var1, 2, mean)
  var1_mean_dt = data.frame(x_axis = rep(x_axis, each = length(y_axis)),
                            y_axis = rep(y_axis, length(x_axis)),
                            z_mean = c(z_var1_mean))
  z_var2_mean = apply(z_var2, 2, mean)
  var2_mean_dt = data.frame(x_axis = rep(x_axis, each = length(y_axis)),
                            y_axis = rep(y_axis, length(x_axis)),
                            z_mean = c(z_var2_mean))
  ## smooth the two-dimensional empirical mean function using the thin plate regression splines
  ## using the 'gam' function in the 'mgcv' package
  ## see '?mgcv::gam' for more details on the configuration of the smoothing process
  var1_mu_fit = as.vector(predict(gam(z_mean ~ s(x_axis, y_axis), data = var1_mean_dt), 
                                  newdata = data.frame(x_axis = rep(x_axis, each = length(y_axis)),
                                                       y_axis = rep(y_axis, length(x_axis)))))
  var2_mu_fit = as.vector(predict(gam(z_mean ~ s(x_axis, y_axis), data = var2_mean_dt), 
                                  newdata = data.frame(x_axis = rep(x_axis, each = length(y_axis)),
                                                       y_axis = rep(y_axis, length(x_axis)))))
  
  
  # fit multilevel univariate FPCA for each single variates
  mod_var1 = mfpca.face_center(Y = z_var1,                          # multilevel univariate (vectorized two-dimensional) functional outcome (var1)
                               id = id_array,                       # subject information for each trial
                               argvals = c(1:M),                    # sampling grids on the two-dimensional functional domain
                               p = 4,                               # degree of b-splines used in FACE algorithm
                               pve = mufpca_pve,                    # percentage of variation explained (pve) threshold for multilevel univariate FPCA 
                                                                    #     (used to determine number of eigencomponents retains at both levels)
                               knots = floor(M/5),                  # number of knots to use for b-splines
                               mu = var1_mu_fit)                    # estimated univariate (vectorized two-dimensional) mean function (var1)
  
  mod_var2 = mfpca.face_center(Y = z_var2,                          # multilevel univariate (vectorized two-dimensional) functional outcome (var2)
                               id = id_array,                       # subject information for each trial
                               argvals = c(1:M),                    # sampling grids on the two-dimensional functional domain
                               p = 4,                               # degree of b-splines used in FACE algorithm
                               pve = mufpca_pve,                    # pve threshold for multilevel univariate FPCA (both levels)
                                                                    #     (used to determine number of eigencomponents retains at both levels)
                               knots = floor(M/5),                  # number of knots to use for b-splines
                               mu = var2_mu_fit)                    # estimated univariate (vectorized two-dimensional) mean function (var2)
  
  # calculate the multivariate eigen components using univariate eigenscores
  # number of univariate eigencomponents at each level for each variate
  m1_lvl1 = mod_var1$npc$level1
  m1_lvl2 = mod_var1$npc$level2
  m2_lvl1 = mod_var2$npc$level1
  m2_lvl2 = mod_var2$npc$level2
  m_sum_lvl1 = m1_lvl1 + m2_lvl1
  m_sum_lvl2 = m1_lvl2 + m2_lvl2
  
  # subject-level multivariate eigencomponents
  ## subject-level univariate eigenscore matrix
  joint_score_lvl1 = cbind(mod_var1$scores$level1[,1:m1_lvl1],
                           mod_var2$scores$level1[,1:m2_lvl1])
  ## subject-level eigenscore covariance matrix
  Z_joint_lvl1 = t(joint_score_lvl1) %*% joint_score_lvl1 / (N-1)
  ## eigen decomposition of subject-level eigenscore covariance matrix
  Z_joint_lvl1.eigen = eigen(Z_joint_lvl1)
  ## subject-level multivariate eigenvalues
  evalues_lvl1_multi_est = Z_joint_lvl1.eigen$values
  ## subject-level multivariate eigenfunctions
  eigen_lvl1_var1_multi_est = mod_var1$efunctions$level1[,1:m1_lvl1] %*%
    Z_joint_lvl1.eigen$vectors[1:m1_lvl1, ]
  eigen_lvl1_var2_multi_est = mod_var2$efunctions$level1[,1:m2_lvl1] %*%
    Z_joint_lvl1.eigen$vectors[c(m1_lvl1 + 1:m2_lvl1), ]
  ## subject-level (shared) multivariate eigenscores
  score_lvl1_multi_est = mod_var1$scores$level1[,1:m1_lvl1] %*%
    Z_joint_lvl1.eigen$vectors[1:m1_lvl1, ] +
    mod_var2$scores$level1[,1:m2_lvl1] %*%
    Z_joint_lvl1.eigen$vectors[c(m1_lvl1 + 1:m2_lvl1), ]
  
  # trial-level multivariate eigencomponents
  ## trial-level univariate eigenscore matrix
  joint_score_lvl2 = cbind(mod_var1$scores$level2[,1:m1_lvl2],
                           mod_var2$scores$level2[,1:m2_lvl2])
  ## trial-level eigenscore covariance matrix
  Z_joint_lvl2 = t(joint_score_lvl2) %*% joint_score_lvl2 / (R_sum-1)
  ## eigen decomposition of trial-level eigenscore covariance matrix
  Z_joint_lvl2.eigen = eigen(Z_joint_lvl2)
  ## trial-level multivariate eigenvalues
  evalues_lvl2_multi_est = Z_joint_lvl2.eigen$values
  ## trial-level multivariate eigenfunctions
  eigen_lvl2_var1_multi_est = mod_var1$efunctions$level2[,1:m1_lvl2] %*%
    Z_joint_lvl2.eigen$vectors[1:m1_lvl2, ]
  eigen_lvl2_var2_multi_est = mod_var2$efunctions$level2[,1:m2_lvl2] %*%
    Z_joint_lvl2.eigen$vectors[c(m1_lvl2 + 1:m2_lvl2), ]
  ## trial-level multivariate eigenscores
  score_lvl2_multi_est = mod_var1$scores$level2[,1:m1_lvl2] %*%
    Z_joint_lvl2.eigen$vectors[1:m1_lvl2, ] +
    mod_var2$scores$level2[,1:m2_lvl2] %*%
    Z_joint_lvl2.eigen$vectors[c(m1_lvl2 + 1:m2_lvl2), ]
  
  # retain a finite number of multivariate eigencomponents based on the pve criteria
  m_multi_pred_lvl1 = min(which(cumsum(evalues_lvl1_multi_est)/sum(evalues_lvl1_multi_est) >= mmfpca_pve))
  m_multi_pred_lvl2 = min(which(cumsum(evalues_lvl2_multi_est)/sum(evalues_lvl2_multi_est) >= mmfpca_pve))
  
  # prediction: reconstruct the multilevel multivariate two-dimensional functional data via 
  # K-L expansion with estimated eigencomponents
  z_multi_pred_var1 = z_multi_pred_var2 = matrix(nrow = M, ncol = R_sum)
  for (i in 1:N){
    if (i == 1){
      start = 0
    }else{
      start = sum(R_list[1:(i-1)])
    }
    for (r in 1:R_list[i]){
      index = start + r
      # prediction using multilevel multivariate MFPCA
      z_multi_pred_var1[index,] = mod_var1$mu +                         # estimated mean function
        score_lvl1_multi_est[i, c(1:m_multi_pred_lvl1)] %*% 
        t(eigen_lvl1_var1_multi_est[,c(1:m_multi_pred_lvl1)]) +         # estimated subject-specific deviation from mean
        score_lvl2_multi_est[index, c(1:m_multi_pred_lvl2)] %*%
        t(eigen_lvl2_var1_multi_est[,c(1:m_multi_pred_lvl2)])           # estimated trial-specific deviation from subject mean
      z_multi_pred_var2[index,] = mod_var2$mu +                         # estimated mean function
        score_lvl1_multi_est[i, c(1:m_multi_pred_lvl1)] %*% 
        t(eigen_lvl1_var2_multi_est[,c(1:m_multi_pred_lvl1)]) +         # estimated subject-specific deviation from mean
        score_lvl2_multi_est[index, c(1:m_multi_pred_lvl2)] %*%
        t(eigen_lvl2_var2_multi_est[,c(1:m_multi_pred_lvl2)])           # estimated trial-specific deviation from subject mean
    }
  }
  
  # store results into list structures
  mu = list()                                                             # mean functions
  mu$var1 = var1_mu_fit
  mu$var2 = var2_mu_fit
  efunctions = list()                                                     # eigenfunctions
  efunctions$lvl1$var1 = eigen_lvl1_var1_multi_est[,1:m_multi_pred_lvl1]
  efunctions$lvl1$var2 = eigen_lvl1_var2_multi_est[,1:m_multi_pred_lvl1]
  efunctions$lvl2$var1 = eigen_lvl2_var1_multi_est[,1:m_multi_pred_lvl2]
  efunctions$lvl2$var2 = eigen_lvl2_var2_multi_est[,1:m_multi_pred_lvl2]
  evalues = list()                                                        # eigenvalues
  evalues$lvl1 = evalues_lvl1_multi_est[1:m_multi_pred_lvl1]
  evalues$lvl2 = evalues_lvl2_multi_est[1:m_multi_pred_lvl2]
  scores = list()                                                         # eigenscores
  scores$lvl1 = score_lvl1_multi_est[,1:m_multi_pred_lvl1]
  scores$lvl2 = score_lvl2_multi_est[,1:m_multi_pred_lvl2]
  z_hat = list()                                                          # trial-specific data reconstruction
  z_hat$var1 = z_multi_pred_var1
  z_hat$var2 = z_multi_pred_var2
  z_data = list()                                                         # original data
  z_data$var1 = z_var1
  z_data$var2 = z_var2
  
  # scale the eigenfunctions, eigenvalues, and eigenscores by sqrt(2) to match the true eigenfunctions
  # if the data are generated using univariate eigenfunctions
  if (! efunctions_multi){
    efunctions$lvl1$var1 = efunctions$lvl1$var1*sqrt(2)
    efunctions$lvl1$var2 = efunctions$lvl1$var2*sqrt(2)
    efunctions$lvl2$var1 = efunctions$lvl2$var1*sqrt(2)
    efunctions$lvl2$var2 = efunctions$lvl2$var2*sqrt(2)
    evalues$lvl1 = evalues$lvl1/2
    evalues$lvl2 = evalues$lvl2/2
    scores$lvl1 = scores$lvl1/sqrt(2)
    scores$lvl2 = scores$lvl2/sqrt(2)
  }
  
  mod_mmfpca = list(mu = mu,
                    efunctions = efunctions,
                    evalues = evalues,
                    scores = scores,
                    z_hat = z_hat,
                    z_data = z_data)
  return(mod_mmfpca)
}






mufpca = function(
    z_matrix,               # multilevel univariate (vectorized two-dimensional) functional data 
                            # (matrix of dimension M*R_sum)
    id_array,               # id information for each trial to identify the subject it belong to (array of length R_sum)
    x_axis,                 # x axis of the two-dimensional functional domain (array of length M_x)
    y_axis,                 # y axis of the two-dimensional functional domain (array of length M_y)
    pve,                    # pve threshold for multilevel univariate FPCA (scalar < 1)
    efunctions_multi = F){  # indicator of whether the data is generated using multivariate eigenfunctions 
                            # (logical, TRUE = using multivariate eigenfunctions, FALSE = using univariate eigenfunctions)

  #########################################################################################
  ## Description: Function that fits multilevel univariate FPCA using fast covariance  
  ##              estimation (FACE) method tailored for two-dimensional functional data 
  ## Definition:  M_x: number of sampling points on the x-axis of the two-dimensional functional domain
  ##              M_y: number of sampling points on the y-axis of the two-dimensional functional domain
  ##              M = M_x * M_y: total number of sampling points on the two-dimensional functional domain
  ##              N: number of subjects
  ##              R_sum: total number of trials
  ## Args:        see above
  ## Returns:     mod_mufpca: a list containing the estimated multilevel univariate 
  ##                          eigencomponents and the corresponding data reconstruction
  #########################################################################################  
  
  N = length(unique(id_array))                         # number of subjects
  R_list = array(dim = N)                              # number of trials for each subject
  for (i in 1:N){
    R_list[i] = sum(as.numeric(id_array == i))
  }
  R_sum = length(id_array)                             # total number of trials
  M = dim(z_matrix)[2]                                 # total number of sampling points
  
  # estimate the smoothed two-dimensional mean function
  ## calculate the empirical mean
  z_mean = apply(z_matrix, 2, mean)
  z_mean_dt = data.frame(x_axis = rep(x_axis, each = length(y_axis)),
                            y_axis = rep(y_axis, length(x_axis)),
                            z_mean = c(z_mean))
  ## smooth the two-dimensional empirical mean function using the thin plate regression splines
  ## using the 'gam' function in the 'mgcv' package
  ## see '?mgcv::gam' for more details on the configuration of the smoothing process
  z_mu_fit = as.vector(predict(gam(z_mean ~ s(x_axis, y_axis), data = z_mean_dt), 
                                  newdata = data.frame(x_axis = rep(x_axis, each = length(y_axis)),
                                                       y_axis = rep(y_axis, length(x_axis)))))
  
  # fit the multilevel univariate FPCA given the smoothed overall mean
  mod_center = mfpca.face_center(Y = z_matrix,                        # multilevel univariate (vectorized two-dimensional) functional outcome
                                 id = id_array,                       # subject information for each trial
                                 argvals = c(1:M),                    # sampling grids on the two-dimensional functional domain
                                 p = 4,                               # degree of b-splines used in FACE algorithm
                                 pve = pve,                           # percentage of variation explained (pve) threshold for multilevel univariate FPCA 
                                                                      #     (used to determine number of eigencomponents retains at both levels)
                                 knots = floor(M/5),                  # number of knots to use for b-splines
                                 mu = z_mu_fit)                       # estimated univariate (vectorized two-dimensional) overall mean function
  
  
  # store results into list structures
  efunctions = list()                                                     # eigenfunctions
  efunctions$lvl1 = mod_center$efunctions$level1
  efunctions$lvl2 = mod_center$efunctions$level2
  evalues = list()                                                        # eigenvalues
  evalues$lvl1 = mod_center$evalues$level1
  evalues$lvl2 = mod_center$evalues$level2
  scores = list()                                                         # eigenscores
  scores$lvl1 = mod_center$scores$level1
  scores$lvl2 = mod_center$scores$level2
  
  # scale the eigenfunctions, eigenvalues, and eigenscores by sqrt(2) to match the true eigenfunctions
  # if the data are generated using multivariate eigenfunctions
  if (efunctions_multi){
    efunctions$lvl1 = efunctions$lvl1/sqrt(2)
    efunctions$lvl2 = efunctions$lvl2/sqrt(2)
    mod_center$evalues$level1 = mod_center$evalues$level1*2
    mod_center$evalues$level2 = mod_center$evalues$level2*2
    mod_center$scores$level1 = mod_center$scores$level1*sqrt(2)
    mod_center$scores$level2 = mod_center$scores$level2*sqrt(2)
  }
  
  # save model components
  mod_mufpca = list(mu = z_mu_fit,                                    # overall mean function
                    efunctions = efunctions, 
                    evalues = evalues, 
                    scores = scores, 
                    z_hat = mod_center$Xhat,                         # original data
                    z_data = z_matrix)                                # trial-specific functional data reconstructed using model components
  
  return(mod_mufpca)
}










mfpca.face_center = function (Y,                    # multilevel univariate functional data (matrix of dimension M*R_sum)       
                              id,                   # id information for each trial to identify the subject it belong to (array of length R_sum)    
                              argvals = NULL,       # sampling grids on the two-dimensional functional domain
                              pve = 0.99,           # pve threshold for multilevel univariate FPCA
                                                    #     (used to determine number of eigencomponents retained at both levels)
                              p = 3,                # degree of freedom of b-splines used in FACE algorithm
                              knots = 35,           # number of knots to use for b-splines
                              mu,                   # overall mean function (array of length R_sum)
                              #####################################################################
                              ## arguments below and their default values are the same as those of
                              ## the 'mfpca.face' function from the 'refund' package
                              ## see '?refund::mfpca.face' for more details
                              #####################################################################
                              visit = NULL, 
                              twoway = FALSE, 
                              weight = "obs", 
                              npc = NULL, 
                              m = 2,
                              silent = TRUE)
{
  
  #########################################################################################
  ## Description: Function that fits multilevel univariate FPCA for one-dimensional functional
  ##              observations using FACE method (Cui et al. 2023) given a pre-estimated mean
  ##              function
  ##              * this function is modified from the 'mfpca.face' function in the 'refund' package
  ##              * the overall mean function needs to be pre-estimated and provided in the input
  ##                (different from the original function)
  ## Definition:  M: total number of sampling points on the two-dimensional functional domain
  ##              N: number of subjects
  ##              R_sum: total number of trials
  ## Args:        see above
  ## Returns:     A list containing
  ##                evalues: estimated subject- and trial-level eigenvalues
  ##                efunctions: estimated subject- and trial-level eigenfunctions
  ##                mu: estimated overall mean function
  ##                scores: estimated subject- and trail-level eigenscores
  ##                Xhat: prediction for trial-level functional data
  ##              see '?refund::mfpca.face' for more details
  ## Note: 1. In line 451-454, the overall mean function is provided by the input 'mu' 
  ##       2. The rest of the function is the same as the 'mfpca.face' function in
  ##          the 'refund' package. See '?refund::mfpca.face' for more details.
  ## Supporting functions: Supporting functions 'face.Cov.mfpca' and 'matrix.multiply.mfpca'
  ##                       utilized in the 'mfpca.face_center' function can be downloaded
  ##                       from https://github.com/refunders/refund/blob/master/R/face.Cov.mfpca.R
  ##                       and manually loaded into the environment by 'source('face.Cov.mfpca.R')'
  #########################################################################################    
  
  
  pspline.setting.mfpca <- function(x, knots = 35, p = 3, 
                                    m = 2, weight = NULL, type = "full", knots.option = "equally-spaced") {
    K = length(knots) - 2 * p - 1
    B = spline.des(knots = knots, x = x, ord = p + 1, outer.ok = TRUE, 
                   sparse = TRUE)$design
    bs = "ps"
    if (knots.option == "quantile") {
      bs = "bs"
    }
    s.object = s(x = x, bs = bs, k = K + p, m = c(p - 1, 
                                                  2), sp = NULL)
    object = smooth.construct(s.object, data = data.frame(x = x), 
                              knots = list(x = knots))
    P = object$S[[1]]
    if (knots.option == "quantile") 
      P = P/max(abs(P)) * 10
    if (is.null(weight)) 
      weight <- rep(1, length(x))
    if (type == "full") {
      Sig = crossprod(matrix.multiply.mfpca(B, weight, 
                                            option = 2), B)
      eSig = eigen(Sig)
      V = eSig$vectors
      E = eSig$values
      if (min(E) <= 1e-07) {
        E <- E + 1e-06
      }
      Sigi_sqrt = matrix.multiply.mfpca(V, 1/sqrt(E)) %*% 
        t(V)
      tUPU = Sigi_sqrt %*% (P %*% Sigi_sqrt)
      Esig = eigen(tUPU, symmetric = TRUE)
      U = Esig$vectors
      s = Esig$values
      s[(K + p - m + 1):(K + p)] = 0
      A = B %*% (Sigi_sqrt %*% U)
    }
    if (type == "simple") {
      A = NULL
      s = NULL
      Sigi_sqrt = NULL
      U = NULL
    }
    List = list(A = A, B = B, s = s, Sigi.sqrt = Sigi_sqrt, 
                U = U, P = P)
    return(List)
  }
  quadWeights.mfpca <- function(argvals, method = "trapezoidal") {
    ret <- switch(method, trapezoidal = {
      D <- length(argvals)
      1/2 * c(argvals[2] - argvals[1], argvals[3:D] - 
                argvals[1:(D - 2)], argvals[D] - argvals[D - 
                                                           1])
    }, midpoint = c(0, diff(argvals)), stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule"))
    return(ret)
  }
  if (silent == FALSE) 
    print("Organize the input")
  stopifnot((!is.null(Y) & !is.null(id)))
  stopifnot(is.matrix(Y))
  if (!is.null(visit)) {
    visit <- as.factor(visit)
  } else {
    visit <- as.factor(ave(id, id, FUN = seq_along))
  }
  id <- as.factor(id)
  df <- data.frame(id = id, visit = visit, Y = I(Y))
  rm(id, visit, Y)
  J <- length(levels(df$visit))
  L <- ncol(df$Y)
  nVisits <- data.frame(table(df$id))
  colnames(nVisits) = c("id", "numVisits")
  ID = sort(unique(df$id))
  I <- length(ID)
  if (is.null(argvals)) 
    argvals <- seq(0, 1, length.out = L)
  if (silent == FALSE) 
    print("Estimate population and visit-specific mean functions")
  # meanY <- colMeans(df$Y, na.rm = TRUE)
  # fit_mu <- gam(meanY ~ s(argvals))
  # mu <- as.vector(predict(fit_mu, newdata = data.frame(argvals = argvals)))
  # rm(meanY, fit_mu)
  mueta = matrix(0, L, J)
  eta = matrix(0, L, J)
  colnames(mueta) <- colnames(eta) <- levels(df$visit)
  Ytilde <- matrix(NA, nrow = nrow(df$Y), ncol = ncol(df$Y))
  if (twoway == TRUE) {
    for (j in 1:J) {
      ind_j <- which(df$visit == levels(df$visit)[j])
      if (length(ind_j) > 1) {
        meanYj <- colMeans(df$Y[ind_j, ], na.rm = TRUE)
      } else {
        meanYj <- df$Y[ind_j, ]
      }
      fit_mueta <- gam(meanYj ~ s(argvals))
      mueta[, j] <- predict(fit_mueta, newdata = data.frame(argvals = argvals))
      eta[, j] <- mueta[, j] - mu
      Ytilde[ind_j, ] <- df$Y[ind_j, ] - matrix(mueta[, 
                                                      j], nrow = length(ind_j), ncol = L, byrow = TRUE)
    }
    rm(meanYj, fit_mueta, ind_j, j)
  } else {
    Ytilde <- df$Y - matrix(mu, nrow = nrow(df$Y), ncol = L, 
                            byrow = TRUE)
  }
  df$Ytilde <- I(Ytilde)
  rm(Ytilde)
  if (silent == FALSE) 
    print("Prepare ingredients for FACE")
  if (length(knots) == 1) {
    if (knots + p >= L) 
      cat("Too many knots!\n")
    stopifnot(knots + p < L)
    K.p <- knots
    knots <- seq(-p, K.p + p, length = K.p + 1 + 2 * p)/K.p
    knots <- knots * (max(argvals) - min(argvals)) + min(argvals)
  }
  if (length(knots) > 1) 
    K.p <- length(knots) - 2 * p - 1
  if (K.p >= L) 
    cat("Too many knots!\n")
  stopifnot(K.p < L)
  c.p <- K.p + p
  List <- pspline.setting.mfpca(argvals, knots, p, m)
  B <- List$B
  Sigi.sqrt <- List$Sigi.sqrt
  s <- List$s
  U <- List$U
  A0 <- Sigi.sqrt %*% U
  G <- crossprod(B)/nrow(B)
  eig_G <- eigen(G, symmetric = T)
  G_half <- eig_G$vectors %*% diag(sqrt(eig_G$values)) %*% 
    t(eig_G$vectors)
  G_invhalf <- eig_G$vectors %*% diag(1/sqrt(eig_G$values)) %*% 
    t(eig_G$vectors)
  Bnew <- as.matrix(B %*% G_invhalf)
  Anew <- G_half %*% A0
  rm(List, Sigi.sqrt, U, G, eig_G, G_half)
  if (silent == FALSE) 
    print("Estimate the total covariance (Kt)")
  Ji <- as.numeric(table(df$id))
  diagD <- rep(Ji, Ji)
  smooth.Gt = face.Cov.mfpca(Y = unclass(df$Ytilde), argvals,                   
                             A0, B, Anew, Bnew, G_invhalf, s)                   
  if (sum(is.na(df$Ytilde)) > 0) {
    df$Ytilde[which(is.na(df$Ytilde))] <- smooth.Gt$Yhat[which(is.na(df$Ytilde))]
  }
  if (weight == "subj") {
    YH <- unclass(df$Ytilde) * sqrt(nrow(df$Ytilde)/(I * 
                                                       diagD))
    smooth.Gt <- face.Cov.mfpca(Y = YH, argvals, A0, B, 
                                Anew, Bnew, G_invhalf, s)
    rm(YH)
  }
  diag_Gt <- colMeans(df$Ytilde^2)
  if (silent == FALSE) 
    print("Estimate principal components of the within covariance (Kw)")
  inx_row_ls <- split(1:nrow(df$Ytilde), f = factor(df$id, 
                                                    levels = unique(df$id)))
  Ysubm <- t(vapply(inx_row_ls, function(x) colSums(df$Ytilde[x, 
                                                              , drop = FALSE], na.rm = TRUE), numeric(L)))
  if (weight == "obs") {
    weights <- sqrt(nrow(df$Ytilde)/(sum(diagD) - nrow(df$Ytilde)))
    YR <- do.call("rbind", lapply(1:I, function(x) {
      weights * sqrt(Ji[x]) * t(t(df$Ytilde[inx_row_ls[[x]], 
                                            , drop = FALSE]) - Ysubm[x, ]/Ji[x])
    }))
  }
  if (weight == "subj") {
    weights <- sqrt(nrow(df$Ytilde)/sum(Ji > 1))
    YR <- do.call("rbind", lapply(1:I, function(x) {
      if (Ji[x] > 1) 
        return((weights/sqrt(Ji[x] - 1)) * t(t(df$Ytilde[inx_row_ls[[x]], 
                                                         , drop = FALSE]) - Ysubm[x, ]/Ji[x]))
    }))
  }
  smooth.Gw <- face.Cov.mfpca(Y = YR, argvals, A0, B, Anew, 
                              Bnew, G_invhalf, s)
  sigma.Gw <- smooth.Gw$evalues
  per <- cumsum(sigma.Gw)/sum(sigma.Gw)
  N.Gw <- ifelse(is.null(npc), min(which(per > pve)), min(npc, 
                                                          length(sigma.Gw)))
  smooth.Gw$efunctions <- smooth.Gw$efunctions[, 1:N.Gw]
  smooth.Gw$evalues <- smooth.Gw$evalues[1:N.Gw]
  rm(Ji, diagD, inx_row_ls, weights, weight, Ysubm, YR, B, 
     Anew, G_invhalf, s, per, N.Gw, sigma.Gw)
  if (silent == FALSE) 
    print("Estimate principal components of the between covariance (Kb)")
  temp = smooth.Gt$decom - smooth.Gw$decom
  Eigen <- eigen(temp, symmetric = TRUE)
  Sigma <- Eigen$values
  d <- Sigma[1:c.p]
  d <- d[d > 0]
  per <- cumsum(d)/sum(d)
  N.Gb <- ifelse(is.null(npc), min(which(per > pve)), min(npc, 
                                                          length(d)))
  smooth.Gb <- list(evalues = Sigma[1:N.Gb], efunctions = Bnew %*% 
                      Eigen$vectors[, 1:N.Gb])
  rm(smooth.Gt, temp, Eigen, Sigma, d, per, N.Gb, Bnew)
  if (silent == FALSE) 
    print("Estimate eigenvalues and eigenfunctions at two levels")
  efunctions <- list(level1 = as.matrix(smooth.Gb$efunctions), 
                     level2 = as.matrix(smooth.Gw$efunctions))
  evalues <- list(level1 = smooth.Gb$evalues, level2 = smooth.Gw$evalues)
  npc <- list(level1 = length(evalues[[1]]), level2 = length(evalues[[2]]))
  names(efunctions) <- names(evalues) <- names(npc) <- c("level1", 
                                                         "level2")
  rm(smooth.Gb, smooth.Gw)
  if (silent == FALSE) 
    print("Estimate the measurement error variance (sigma^2)")
  cov.hat <- lapply(c("level1", "level2"), function(x) colSums(t(efunctions[[x]]^2) * 
                                                                 evalues[[x]]))
  T.len <- argvals[L] - argvals[1]
  T1.min <- min(which(argvals >= argvals[1] + 0.25 * T.len))
  T1.max <- max(which(argvals <= argvals[L] - 0.25 * T.len))
  DIAG <- (diag_Gt - cov.hat[[1]] - cov.hat[[2]])[T1.min:T1.max]
  w2 <- quadWeights.mfpca(argvals[T1.min:T1.max], method = "trapezoidal")
  sigma2 <- max(weighted.mean(DIAG, w = w2, na.rm = TRUE), 
                0)
  rm(cov.hat, T.len, T1.min, T1.max, DIAG, w2)
  if (silent == FALSE) 
    print("Estimate principal component scores")
  Xhat <- Xhat.subject <- matrix(0, nrow(df$Y), L)
  phi1 <- efunctions[[1]]
  phi2 <- efunctions[[2]]
  score1 <- matrix(0, I, npc[[1]])
  score2 <- matrix(0, nrow(df$Y), npc[[2]])
  unVisits <- unique(nVisits$numVisits)
  if (length(unVisits) < I) {
    for (j in 1:length(unVisits)) {
      Jm <- unVisits[j]
      if (sigma2 < 1e-04) {
        A <- Jm * (t(phi1) %*% phi1)
        B <- matrix(rep(t(phi1) %*% phi2, Jm), nrow = npc[[1]])
        temp <- ginv(t(phi2) %*% phi2)
      } else {
        if (length(evalues[[1]]) == 1) {
          A <- Jm * (t(phi1) %*% phi1)/sigma2 + 1/evalues[[1]]
        } else {
          A <- Jm * (t(phi1) %*% phi1)/sigma2 + diag(1/evalues[[1]])
        }
        B = matrix(rep(t(phi1) %*% phi2/sigma2, Jm), 
                   nrow = npc[[1]])
        if (length(evalues[[2]]) == 1) {
          temp = ginv(t(phi2) %*% phi2/sigma2 + 1/evalues[[2]])
        } else {
          temp = ginv(t(phi2) %*% phi2/sigma2 + diag(1/evalues[[2]]))
        }
      }
      C <- t(B)
      invD <- kronecker(diag(1, Jm), temp)
      MatE <- ginv(A - B %*% invD %*% C)
      MatF <- -invD %*% C %*% MatE
      MatG <- -MatE %*% B %*% invD
      MatH <- invD - invD %*% C %*% MatG
      Mat1 <- cbind(MatE, MatG)
      Mat2 <- cbind(MatF, MatH)
      ind.Jm <- nVisits$id[which(nVisits$numVisits == 
                                   Jm)]
      YJm <- matrix(df$Ytilde[which(df$id %in% ind.Jm), 
      ], ncol = L)
      int1 <- rowsum(df$Ytilde[which(df$id %in% ind.Jm), 
      ] %*% phi1, rep(1:length(ind.Jm), each = Jm))
      int2 <- t(matrix(t(df$Ytilde[which(df$id %in% ind.Jm), 
      ] %*% phi2), nrow = npc[[2]] * Jm))
      int <- cbind(int1, int2)
      if (sigma2 >= 1e-04) {
        int <- int/sigma2
      }
      score1[which(nVisits$id %in% ind.Jm), ] <- int %*% 
        t(Mat1)
      score2[which(df$id %in% ind.Jm), ] <- t(matrix(Mat2 %*% 
                                                       t(int), nrow = npc[[2]]))
      temp <- score1[which(nVisits$id %in% ind.Jm), ] %*% 
        t(phi1)
      Xhat.subject[which(df$id %in% ind.Jm), ] <- temp[rep(1:length(ind.Jm), 
                                                           each = Jm), ]
      Xhat[which(df$id %in% ind.Jm), ] <- Xhat.subject[which(df$id %in% 
                                                               ind.Jm), ] + score2[which(df$id %in% ind.Jm), 
                                                               ] %*% t(phi2)
    }
    for (g in 1:length(levels(df$visit))) {
      ind.visit <- which(df$visit == levels(df$visit)[g])
      Xhat.subject[ind.visit, ] <- t(t(Xhat.subject[ind.visit, 
      ]) + mu + eta[, levels(df$visit)[g]])
      Xhat[ind.visit, ] <- t(t(Xhat[ind.visit, ]) + mu + 
                               eta[, levels(df$visit)[g]])
    }
    rm(YJm, g, ind.visit, ind.Jm)
  } else {
    for (m in 1:I) {
      Jm <- nVisits[m, 2]
      if (sigma2 < 1e-04) {
        A <- Jm * (t(phi1) %*% phi1)
        B <- matrix(rep(t(phi1) %*% phi2, Jm), nrow = npc[[1]])
        temp <- ginv(t(phi2) %*% phi2)
      } else {
        if (length(evalues[[1]]) == 1) {
          A <- Jm * (t(phi1) %*% phi1)/sigma2 + 1/evalues[[1]]
        } else {
          A <- Jm * (t(phi1) %*% phi1)/sigma2 + diag(1/evalues[[1]])
        }
        B = matrix(rep(t(phi1) %*% phi2/sigma2, Jm), 
                   nrow = npc[[1]])
        if (length(evalues[[2]]) == 1) {
          temp = ginv(t(phi2) %*% phi2/sigma2 + 1/evalues[[2]])
        }  else {
          temp = ginv(t(phi2) %*% phi2/sigma2 + diag(1/evalues[[2]]))
        }
      }
      C <- t(B)
      invD <- kronecker(diag(1, Jm), temp)
      MatE <- ginv(A - B %*% invD %*% C)
      MatF <- -invD %*% C %*% MatE
      MatG <- -MatE %*% B %*% invD
      MatH <- invD - invD %*% C %*% MatG
      Mat1 <- cbind(MatE, MatG)
      Mat2 <- cbind(MatF, MatH)
      int1 <- colSums(matrix(df$Ytilde[df$id == ID[m], 
      ], ncol = L) %*% phi1)
      int2 <- matrix(df$Ytilde[df$id == ID[m], ], ncol = L) %*% 
        phi2
      if (sigma2 < 1e-04) {
        int <- c(int1, as.vector(t(int2)))
      } else {
        int <- c(int1, as.vector(t(int2)))/sigma2
      }
      score1[m, ] <- Mat1 %*% int
      score2[which(df$id == ID[m]), ] <- matrix(Mat2 %*% 
                                                  int, ncol = npc[[2]], byrow = TRUE)
      for (j in which(df$id == ID[m])) {
        Xhat.subject[j, ] <- as.matrix(mu) + eta[, df$visit[j]] + 
          as.vector(phi1 %*% score1[m, ])
        Xhat[j, ] <- Xhat.subject[j, ] + as.vector(phi2 %*% 
                                                     score2[j, ])
      }
    }
  }
  scores <- list(level1 = score1, level2 = score2)
  rm(A, B, C, int, int1, int2, invD, Mat1, Mat2, MatE, MatF, 
     MatG, MatH, temp, j, Jm, unVisits, phi1, phi2, score1, 
     score2)
  if (silent == FALSE) 
    print("Organize the results")
  res <- list(Xhat = Xhat, Xhat.subject = Xhat.subject, Y = df$Y, 
              mu = mu, eta = eta, scores = scores, efunctions = efunctions, 
              evalues = evalues, npc = npc, sigma2 = sigma2)
  rm(df, efunctions, eta, evalues, mueta, nVisits, npc, scores, 
     Xhat, Xhat.subject, argvals, diag_Gt, I, ID, J, mu, 
     pve, L, sigma2)
  return(res)
}



