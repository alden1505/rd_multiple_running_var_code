#' Subfunction that creates basic MRD estimates
#' 
#' basic_est is a subfunction that does the basic MRD estimation
#' 
#' @import locfit
#' @import mgcv
#' @import Hmisc
#' @import ggplot2
#' @param Y Y is a vector for the dependent variable.
#' @param X1 X1 is a vector for the first running variable.
#' @param X2 X2 is a vector for the second running variable.
#' @param W W is a vector containing dummy variables for whether the observation is treated. This defaults to NULL, and it is assumed to be a sharp RD where a unit is treated if and only if both running variables exceed their threshold. If it is in fact a fuzzy RD (so the fuzzy option is specified to be TRUE), then W must be specified.
#' @param F1.range This is the points at which the user wants to compute the conditional average treatment effect at the boundary where X2=0, and X1 is positive. The default for this option is NULL, in which case the conditional average treatment effect is computed from X1=0 to max(X1), at equal intervals.
#' @param F2.range This is the points at which the user wants to compute the conditional average treatment effect at the boundary where X1=0, and X2 is positive. The default for this option is NULL, in which case the conditional average treatment effect is computed from X2=0 to max(X2), at equal intervals.
#' @param R1 The input for this option should be a vector containing indices for observations in treated region. Typically, this can be left blank, in which case R1 defaults to NULL.
#' @param R2 The input for this option should be a vector containing indices for observations in untreated region. Typically, this can be left blank, in which case R2 defaults to NULL.
#' @param K This is a parameter controlling the penalization for thin plate regression splines. Typically this is left empty, in which case K defaults to -1.
#' @param save_bayes_se The input should be a logical value (TRUE or FALSE) indicating whether Bayesian standard errors should be saved.
#' @param fuzzy The input should be a logical value (TRUE or FALSE) indicating whether this is a fuzzy RD. If this is not specified, the option defaults to FALSE (so the code assumes that this is a sharp RD).
#' @param reverse The input should be a logical value (TRUE or FALSE) indicating whether the running variables are reversed (i.e. negative values are treated). If not specified, the option defaults to FALSE.
#' @param sp This is a parameter for thin plate regression spline fitting, and can typically be left empty (in which case it defaults to NULL).
#' @param save_fits The input should be a logical value indicating whether the fitted thin plate regression splines should be saved. If the option is not specified, it defaults to TRUE.
#' @param save_density_est The input should be a logical value indicating whether kernel density estimates should be saved. If the option is not specified, it defaults to TRUE.
#' @param pred_dfs The input should be a dataframe with the values of the running variables that conditional average treatment effects should be calculated at. If left empty, it defaults to NULL.
#' @param undersmooth The input should be a logical value indicating whether thin plate regression splines should be undersmoothed, in order to reduce bias (for inference purposes).
#' @param save_sp The input should be a logical value indicating whether the parameter for thin plate regression spline fitting should be saved. If not specified, it defaults to TRUE.
#' @keywords basic_mrd_est
#' @export

basic_est = function(Y, X1, X2, W = NULL, F1.range = NULL, F2.range = NULL,
                     R1 = NULL, R2 = NULL, K = -1, save_bayes_se = TRUE,
                     fuzzy = FALSE, reverse = FALSE, sp = NULL,
                     save_fits = TRUE, save_density_est = TRUE,
                     pred_dfs = NULL, undersmooth, save_sp = TRUE){
  ### Set F1 and F2 if they were not specified
  if (is.null(F1.range)){
    F1.range = seq(0, max(X2), length = 100)
  }
  if (is.null(F2.range)){
    F2.range = seq(0, max(X1), length = 100)
  }
  ### Find indices of observations for the 2 regions
  if (is.null(R1)){
    R1 = which(X1>=0 & X2>=0)
  }
  if (is.null(R2)){
    R2 = which(X1<0 | X2<0)
  }
  regions = list(R1,R2)
  if (is.null(pred_dfs)){
    pred_dfs = list(data.frame(X1=0, X2=F1.range),
                    data.frame(X2=0, X1=F2.range))
  }
  # Do density estimation
  kernel_est = list(preds = list(), mass = list())
  kernel_density_est = locfit(~ X1 + X2)
  for (f in 1:2){
    kernel_est$preds[[f]] = predict(kernel_density_est, pred_dfs[[f]])
    kernel_est$mass[[f]] = kernel_est$preds[[f]]/sum(kernel_est$preds[[f]])
  }
  kernel_est$mass[[3]] = c(kernel_est$preds[[1]],kernel_est$preds[[2]])/
    sum(c(kernel_est$preds[[1]],kernel_est$preds[[2]]))
  
  # Total number of effects (e.g. reduced form, first stage, LATE for fuzzy RD)
  if (fuzzy){
    E = 3
    fits = list(reduced = list(), stage1 = list())
  }else{
    E = 1
    fits = list(reduced = list())
  }
  if (undersmooth){
    fits_us = fits
  }
  
  # Fit splines
  for (r in 1:2){
    if (is.null(sp)){
      fits$reduced[[r]] = gam(Y ~ s(X1, X2, k = K), subset = regions[[r]])
      if (undersmooth){
        fits_us$reduced[[r]] = gam(Y ~ s(X1, X2, k = K), subset = regions[[r]],
                                   sp = 0.5*fits$reduced[[r]]$sp)
      }
      if (fuzzy){
        fits$stage1[[r]] = gam(W ~ s(X1, X2, k = K), subset = regions[[r]])
        if (undersmooth){
          fits_us$stage1[[r]] = gam(W ~ s(X1, X2, k = K), subset = regions[[r]],
                                    sp = 0.5*fits$stage1[[r]]$sp)
        }
      }
    }else{
      fits$reduced[[r]] = gam(Y ~ s(X1, X2, k = K), subset = regions[[r]], sp = sp$reduced[r])
      if (undersmooth){
        fits_us$reduced[[r]] = gam(Y ~ s(X1, X2, k = K), subset = regions[[r]],
                                   sp = 0.5*sp$reduced[r])
      }
      if (fuzzy){
        fits$stage1[[r]] = gam(W ~ s(X1, X2, k = K), subset = regions[[r]], sp = sp$stage1[r])
        if (undersmooth){
          fits_us$stage1[[r]] = gam(W ~ s(X1, X2, k = K), subset = regions[[r]],
                                    sp = 0.5*sp$stage1[r])
        }
      }
    }
  }
  # Obtain predictions
  if (fuzzy){
    preds = list(reduced = list(F1 = list(), F2 = list()),
                 stage1 = list(F1 = list(), F2 = list()))
  }else{
    preds = list(reduced = list(F1 = list(), F2 = list()))
  }
  preds_us = preds
  for (e in 1:max(E-1,1)){
    for (f in 1:2){
      for (r in 1:2){
        if (save_bayes_se){
          preds[[e]][[f]][[r]] = predict(fits[[e]][[r]], pred_dfs[[f]], se.fit = TRUE)
          if (undersmooth){
            preds_us[[e]][[f]][[r]] = predict(fits_us[[e]][[r]], pred_dfs[[f]], se.fit = TRUE)
          }
        }else{
          preds[[e]][[f]][[r]] = predict(fits[[e]][[r]], pred_dfs[[f]])
          if (undersmooth){
            preds_us[[e]][[f]][[r]] = predict(fits_us[[e]][[r]], pred_dfs[[f]])
          }
        }
      }
    }
  }
  # Treatment effects
  tau = list(reduced = list(pw = list(F1 = list(), F2 = list()),
                            agg = list(F1 = list(), F2 = list(), all = list())))
  if (fuzzy){
    tau$stage1 = list(pw = list(F1 = list(), F2 = list()),
                      agg = list(F1 = list(), F2 = list(), all = list()))
    tau$late = list(pw = list(F1 = list(), F2 = list()),
                    agg = list(F1 = list(), F2 = list(), all = list()))
  }
  if (undersmooth){
    tau_us = tau
  }
  for (f in 1:2){
    for (e in 1:max(E-1,1)){
      if (save_bayes_se){
        tau[[e]]$pw[[f]] = preds[[e]][[f]][[1]]$fit - preds[[e]][[f]][[2]]$fit
        if (undersmooth){
          tau_us[[e]]$pw[[f]] = preds_us[[e]][[f]][[1]]$fit - preds_us[[e]][[f]][[2]]$fit
        }
      }else{
        tau[[e]]$pw[[f]] = preds[[e]][[f]][[1]] - preds[[e]][[f]][[2]]
        if (undersmooth){
          tau_us[[e]]$pw[[f]] = preds_us[[e]][[f]][[1]] - preds_us[[e]][[f]][[2]]
        }
      }
      if (reverse){
        tau[[e]]$pw[[f]] = -tau[[e]]$pw[[f]]
        if (undersmooth){
          tau_us[[e]]$pw[[f]] = -tau_us[[e]]$pw[[f]]
        }
      }
    }
    if (fuzzy){
      tau$late$pw[[f]] = tau$reduced$pw[[f]]/tau$stage1$pw[[f]]
      if (undersmooth){
        tau_us$late$pw[[f]] = tau_us$reduced$pw[[f]]/tau_us$stage1$pw[[f]]
      }
    }
    for (e in 1:E){
      tau[[e]]$agg[[f]] = sum(tau[[e]]$pw[[f]]*kernel_est$mass[[f]])
      if (undersmooth){
        tau_us[[e]]$agg[[f]] = sum(tau_us[[e]]$pw[[f]]*kernel_est$mass[[f]])
      }
    }
  }
  for (e in 1:E){
    tau[[e]]$agg[[3]] = sum(c(tau[[e]]$pw[[1]],tau[[e]]$pw[[2]])*
                              kernel_est$mass[[3]])
    if (undersmooth){
      tau_us[[e]]$agg[[3]] = sum(c(tau_us[[e]]$pw[[1]],tau_us[[e]]$pw[[2]])*
                                   kernel_est$mass[[3]])
    }
  }
  # Save Bayesian standard errors if desired
  if (save_bayes_se){
    tau_se = list(reduced = list(pw = list(bayes = list()),
                                 agg=list(delta = list(F1=list(),F2=list()))))
    if (fuzzy){
      tau_se$stage1 = list(pw = list(bayes = list()),
                           agg=list(delta = list(F1=list(),F2=list())))
    }
    if (undersmooth){
      tau_se_us = tau_se
    }
    for (e in 1:max(E-1,1)){
      tau_se[[e]]$pw$bayes$cor = list(F1 = list(), F2 = list())
      tau_se[[e]]$pw$bayes$ind = list(F1 = list(), F2 = list())
      if (undersmooth){
        tau_se_us$pw = tau_se
      }
      for (f in 1:2){
        tau_se[[e]]$pw$bayes$cor[[f]] =
          sqrt(preds[[e]][[f]][[1]]$se.fit^2 +
                 preds[[e]][[f]][[2]]$se.fit^2 +
                 2*preds[[e]][[f]][[1]]$se.fit*preds[[e]][[f]][[2]]$se.fit)
        tau_se[[e]]$pw$bayes$ind[[f]] =
          sqrt(preds[[e]][[f]][[1]]$se.fit^2 +
                 preds[[e]][[f]][[2]]$se.fit^2)
        if (undersmooth){
          tau_se_us[[e]]$pw$bayes$cor[[f]] =
            sqrt(preds_us[[e]][[f]][[1]]$se.fit^2 +
                   preds_us[[e]][[f]][[2]]$se.fit^2 +
                   2*preds_us[[e]][[f]][[1]]$se.fit*preds_us[[e]][[f]][[2]]$se.fit)
          tau_se_us[[e]]$pw$bayes$ind[[f]] =
            sqrt(preds_us[[e]][[f]][[1]]$se.fit^2 +
                   preds_us[[e]][[f]][[2]]$se.fit^2)
        }
      }
      
    }
  }
  results = list(tau = tau)
  if (undersmooth){
    results$tau_us = tau_us
  }
  if (save_bayes_se){
    results$tau_se = tau_se
    if (undersmooth){
      results$tau_se_us = tau_se_us
    }
  }
  if (save_fits){
    results$fits = fits
    if (undersmooth){
      results$fits_us = fits_us
    }
  }
  if (save_density_est){
    results$density_est = kernel_est$mass
    results$density_preds = kernel_est$preds
  }
  if (save_sp){
    results$sp = list(reduced = c(fits$reduced[[1]]$sp,fits$reduced[[2]]$sp))
    if (fuzzy){
      results$sp$stage1 = c(fits$stage1[[1]]$sp,fits$stage1[[2]]$sp)
    }
  }
  return(results)
}

#' Main function for MRD Estimation
#' 
#' mrdd.fn is the main function for MRD estimation
#' 
#' @import locfit
#' @import mgcv
#' @import Hmisc
#' @import ggplot2
#' @param Y Y is a vector for the dependent variable.
#' @param X1 X1 is a vector for the first running variable.
#' @param X2 X2 is a vector for the second running variable.
#' @param W W is a vector containing dummy variables for whether the observation is treated. This defaults to NULL, and it is assumed to be a sharp RD where a unit is treated if and only if both running variables exceed their threshold. If it is in fact a fuzzy RD (so the fuzzy option is specified to be TRUE), then W must be specified.
#' @param F1.range This is the points at which the user wants to compute the conditional average treatment effect at the boundary where X2=0, and X1 is positive. The default for this option is NULL, in which case the conditional average treatment effect is computed from X1=0 to max(X1), at equal intervals.
#' @param F2.range This is the points at which the user wants to compute the conditional average treatment effect at the boundary where X1=0, and X2 is positive. The default for this option is NULL, in which case the conditional average treatment effect is computed from X2=0 to max(X2), at equal intervals.
#' @param R1 The input for this option should be a vector containing indices for observations in treated region. Typically, this can be left blank, in which case R1 defaults to NULL.
#' @param R2 The input for this option should be a vector containing indices for observations in untreated region. Typically, this can be left blank, in which case R2 defaults to NULL.
#' @param zscore The input should be the z-score that should be used for statistical significance. If not specified, it defaults to qnorm(0.975) i.e. 5 percent significance level.
#' @param K This is a parameter controlling the penalization for thin plate regression splines. Typically this is left empty, in which case K defaults to -1.
#' @param bs The input should be a logical value indicating whether to use bootstrap for standard errors. If not specified, this option defaults to FALSE.
#' @param bs_reps The input should be the number of desired bootstrap iterations. If not specified, this is option defaults to 1000.
#' @param bs_seeds The input should be a logical value indicating whether the random seed for the rth bootstrap iteration should be set to r. If not specified, this option defaults to TRUE.
#' @param save_bayes_se The input should be a logical value (TRUE or FALSE) indicating whether Bayesian standard errors should be saved.
#' @param save_density_est The input should be a logical value indicating whether kernel density estimates should be saved. If the option is not specified, it defaults to TRUE.
#' @param save_fits The input should be a logical value indicating whether the fitted thin plate regression splines should be saved. If the option is not specified, it defaults to TRUE.
#' @param fuzzy The input should be a logical value (TRUE or FALSE) indicating whether this is a fuzzy RD. If this is not specified, the option defaults to FALSE (so the code assumes that this is a sharp RD).
#' @param reverse The input should be a logical value (TRUE or FALSE) indicating whether the running variables are reversed (i.e. negative values are treated). If not specified, the option defaults to FALSE.
#' @param undersmooth The input should be a logical value indicating whether thin plate regression splines should be undersmoothed, in order to reduce bias (for inference purposes).
#' @param pred_dfs The input should be a dataframe with the values of the running variables that conditional average treatment effects should be calculated at. If left empty, it defaults to NULL.

mrdd.fn = function(Y, X1, X2, W = NULL, F1.range = NULL, F2.range = NULL,
                   R1 = NULL, R2 = NULL,
                   zscore=qnorm(0.975), K = -1,
                   bs = FALSE, bs_reps = 1000, bs_seeds = TRUE,
                   save_bayes_se = TRUE, save_density_est = TRUE, save_fits = TRUE,
                   fuzzy = FALSE, reverse = FALSE,
                   undersmooth = TRUE,
                   pred_dfs = NULL){
  # Main estimates
  main_est = basic_est(Y = Y, X1 = X1, X2 = X2, W = W, K=K,
                       F1.range = F1.range, F2.range = F2.range, R1 = R1, R2 = R2,
                       save_bayes_se=save_bayes_se, save_density_est = save_density_est,
                       fuzzy=fuzzy, reverse = reverse,
                       undersmooth = undersmooth, save_sp = TRUE,
                       pred_dfs = pred_dfs)
  # Total number of effects (e.g. reduced form, first stage, LATE for fuzzy RD)
  if (fuzzy){
    E = 3
  }else{
    E = 1
  }
  if (is.null(pred_dfs)){
    pred_dfs = list(data.frame(X1=0, X2=F1.range),
                    data.frame(X2=0, X1=F2.range))
  }
  # Create matrices to store bootstrap results
  if (bs | fuzzy){
    bs_est = list(reduced = list(pw = list(F1 = list(), F2 = list()),
                                 agg = list(F1 = list(), F2 = list(), all = list())))
    if (fuzzy){
      bs_est$stage1 = list(pw = list(F1 = list(), F2 = list()),
                           agg = list(F1 = list(), F2 = list(), all = list()))
      bs_est$late = list(pw = list(F1 = list(), F2 = list()),
                         agg = list(F1 = list(), F2 = list(), all = list()))
    }
    for (e in 1:E){
      for (f in 1:2){
        bs_est[[e]]$pw[[f]] =
          matrix(nrow = length(main_est$tau$reduced$pw[[f]]),ncol = bs_reps)
        bs_est[[e]]$agg[[f]] = rep(NA, bs_reps)
      }
      bs_est[[e]]$agg[[3]] = rep(NA, bs_reps)
    }
    if (save_density_est){
      bs_est$density_est = list(F1 = matrix(nrow = length(F1.range), ncol = bs_reps),
                                F2 = matrix(nrow = length(F2.range), ncol = bs_reps),
                                all = matrix(nrow = length(F1.range) + length(F2.range), ncol = bs_reps))
      bs_est$density_preds = list(F1 = matrix(nrow = length(F1.range), ncol = bs_reps),
                                  F2 = matrix(nrow = length(F2.range), ncol = bs_reps),
                                  all = matrix(nrow = length(F1.range) + length(F2.range), ncol = bs_reps))
    }
  }else{
    bs_est = list(density_est = list(F1 = matrix(nrow = length(F1.range), ncol = bs_reps),
                                     F2 = matrix(nrow = length(F2.range), ncol = bs_reps),
                                     all = matrix(nrow = length(F1.range) + length(F2.range), ncol = bs_reps)),
                  density_preds = list(F1 = matrix(nrow = length(F1.range), ncol = bs_reps),
                                       F2 = matrix(nrow = length(F2.range), ncol = bs_reps),
                                       all = matrix(nrow = length(F1.range) + length(F2.range), ncol = bs_reps)))
  }
  if (undersmooth){
    bs_est_us = bs_est
  }
  
  for (b in 1:bs_reps){
    if (!(b%%100)){
      print(b)
    }
    if (bs_seeds){
      set.seed(b)
    }
    bs_indices = sample.int(n = length(Y), size = length(Y), replace = TRUE)
    bs_X1 = X1[bs_indices]
    bs_X2 = X2[bs_indices]
    if (bs){
      bs_Y = Y[bs_indices]
      current_bs_est = basic_est(Y = bs_Y, X1 = bs_X1, X2 = bs_X2,
                                 F1.range = F1.range, F2.range = F2.range,
                                 save_bayes_se = FALSE, save_density_est = save_density_est,
                                 fuzzy=fuzzy, undersmooth = undersmooth, sp = main_est$sp,
                                 pred_dfs = pred_dfs)
      if (fuzzy){
        bs_W = W[bs_indices]
        current_bs_est = basic_est(Y = bs_Y, X1 = bs_X1, X2 = bs_X2, W = bs_W,
                                   F1.range = F1.range, F2.range=F2.range,
                                   save_bayes_se = FALSE, save_density_est = save_density_est,
                                   fuzzy = fuzzy, undersmooth = undersmooth, sp = main_est$sp,
                                   pred_dfs = pred_dfs)
      }
      for (e in 1:E){
        for (f in 1:2){
          bs_est[[e]]$pw[[f]][,b] = current_bs_est$tau[[e]]$pw[[f]]
          bs_est[[e]]$agg[[f]][b] = current_bs_est$tau[[e]]$agg[[f]]
          if (undersmooth){
            bs_est_us[[e]]$pw[[f]][,b] = current_bs_est$tau_us[[e]]$pw[[f]]
            bs_est_us[[e]]$agg[[f]][b] = current_bs_est$tau_us[[e]]$agg[[f]]
          }
        }
        bs_est[[e]]$agg[[3]][b] = current_bs_est$tau[[e]]$agg[[3]]
        if (undersmooth){
          bs_est_us[[e]]$agg[[3]][b] = current_bs_est$tau_us[[e]]$agg[[3]]
        }
      }
    }
    # Do density estimation
    if (is.null(pred_dfs)){
      bs_pred_dfs = list(data.frame(bs_X1=0, bs_X2=F1.range),
                         data.frame(bs_X2=0, bs_X1=F2.range))
    }else{
      bs_pred_dfs = list(data.frame(bs_X1 = pred_dfs[[1]]$X1, bs_X2 = pred_dfs[[1]]$X2),
                         data.frame(bs_X1 = pred_dfs[[2]]$X1, bs_X2 = pred_dfs[[2]]$X2))
    }
    kernel_est = list(preds = list(), mass = list())
    kernel_density_est = locfit(~ bs_X1 + bs_X2)
    for (f in 1:2){
      kernel_est$preds[[f]] = predict(kernel_density_est, bs_pred_dfs[[f]])
      kernel_est$mass[[f]] = kernel_est$preds[[f]]/sum(kernel_est$preds[[f]])
    }
    kernel_est$mass[[3]] = c(kernel_est$preds[[1]],kernel_est$preds[[2]])/
      sum(c(kernel_est$preds[[1]],kernel_est$preds[[2]]))
    
    if (save_density_est){
      for (f in 1:3){
        bs_est$density_est[[f]][,b] = kernel_est$mass[[f]]
        if (f != 3){
          bs_est$density_preds[[f]][,b] = kernel_est$preds[[f]]
        }
      }
    }
  }
  tau_se = main_est$tau_se
  if (undersmooth){
    tau_se_us = main_est$tau_se_us
  }
  if (fuzzy){
    tau_se$late = list(pw = list(), agg = list())
    if (undersmooth){
      tau_se_us$late = list(pw = list(), agg = list())
    }
  }
  # Compute Bayesian SEs for aggregate effects
  var_mat = list(reduced = list(pw = list(F1 = list(), F2 = list())),
                 density_est = list(pw = list(F1 = list(), F2 = list())))
  if (fuzzy){
    var_mat$stage1 = list(pw = list(F1 = list(), F2 = list()))
  }
  if (undersmooth){
    var_mat_us = var_mat
  }
  for (e in 1:max(E-1,1)){
    for (f in 1:2){
      # covariance matrix of treatment effect point estimates
      var_mat[[e]][[f]] = 
        predict(main_est$fits[[e]][[1]], pred_dfs[[f]], type = "lpmatrix") %*% 
        main_est$fits[[e]][[1]]$Vp %*%
        t(predict(main_est$fits[[e]][[1]], pred_dfs[[f]], type = "lpmatrix")) + 
        predict(main_est$fits[[e]][[2]], pred_dfs[[f]], type = "lpmatrix") %*% 
        main_est$fits[[e]][[2]]$Vp %*%
        t(predict(main_est$fits[[e]][[2]], pred_dfs[[f]], type = "lpmatrix"))
      # covariance matrix for density estimate
      var_mat$density_est[[f]] = var(t(bs_est$density_preds[[f]]))/length(Y)
      Deriv_wrt_f = as.numeric( (1/(sum(main_est$density_preds[[f]]))^2) * 
                                  ( sum(main_est$density_preds[[f]]) * main_est$tau[[e]]$pw[[f]] - 
                                      sum(main_est$density_preds[[f]] * main_est$tau[[e]]$pw[[f]]) ) )
      Deriv_wrt_tau = main_est$density_est[[f]]
      tau_se[[e]]$agg$delta[[f]] = sqrt(t(Deriv_wrt_tau) %*% var_mat[[e]][[f]] %*% Deriv_wrt_tau + 
                                          t(Deriv_wrt_f) %*% var_mat$density_est[[f]] %*% Deriv_wrt_f)
      if (undersmooth){
        # covariance matrix of treatment effect point estimates
        var_mat_us[[e]][[f]] = 
          predict(main_est$fits_us[[e]][[1]], pred_dfs[[f]], type = "lpmatrix") %*% 
          main_est$fits_us[[e]][[1]]$Vp %*%
          t(predict(main_est$fits_us[[e]][[1]], pred_dfs[[f]], type = "lpmatrix")) + 
          predict(main_est$fits_us[[e]][[2]], pred_dfs[[f]], type = "lpmatrix") %*% 
          main_est$fits_us[[e]][[2]]$Vp %*%
          t(predict(main_est$fits_us[[e]][[2]], pred_dfs[[f]], type = "lpmatrix"))
        # covariance matrix for density estimate
        var_mat_us$density_est[[f]] = var(t(bs_est$density_preds[[f]]))/length(Y)
        Deriv_wrt_f = as.numeric( (1/(sum(main_est$density_preds[[f]]))^2) * 
                                    ( sum(main_est$density_preds[[f]]) * main_est$tau_us[[e]]$pw[[f]] - 
                                        sum(main_est$density_preds[[f]] * main_est$tau_us[[e]]$pw[[f]]) ) )
        Deriv_wrt_tau = main_est$density_est[[f]]
        tau_se_us[[e]]$agg$delta[[f]] = sqrt(t(Deriv_wrt_tau) %*% var_mat_us[[e]][[f]] %*% Deriv_wrt_tau + 
                                               t(Deriv_wrt_f) %*% var_mat_us$density_est[[f]] %*% Deriv_wrt_f)
      }
    }
  }
  
  if (bs | fuzzy){
    for (e in 1:E){
      tau_se[[e]]$pw$bs = list(F1 = list(), F2 = list())
      tau_se[[e]]$agg$bs = list(F1 = list(), F2 = list(), all = list())
      if (undersmooth){
        tau_se_us[[e]]$pw$bs = list(F1 = list(), F2 = list())
        tau_se_us[[e]]$agg$bs = list(F1 = list(), F2 = list(), all = list())
      }
      for (f in 1:2){
        tau_se[[e]]$pw$bs[[f]] = apply(bs_est[[e]]$pw[[f]], 1, sd)
        tau_se[[e]]$agg$bs[[f]] = sd(bs_est[[e]]$agg[[f]])
        if (undersmooth){
          tau_se_us[[e]]$pw$bs[[f]] = apply(bs_est_us[[e]]$pw[[f]], 1, sd)
          tau_se_us[[e]]$agg$bs[[f]] = sd(bs_est_us[[e]]$agg[[f]])
        }
      }
      tau_se[[e]]$agg$bs[[3]] = sd(bs_est[[e]]$agg[[3]])
      if (undersmooth){
        tau_se_us[[e]]$agg$bs[[3]] = sd(bs_est_us[[e]]$agg[[3]])
      }
    }
    results = list(tau = main_est$tau, tau_se = tau_se, bs_est = bs_est)
  }else{
    results = list(tau = main_est$tau, tau_se = tau_se)
  }
  if (undersmooth){
    results$tau_us = main_est$tau_us
    results$tau_se_us = tau_se_us
    if (bs | fuzzy){
      results$bs_est_us = bs_est_us
    }
  }
  if (save_fits){
    results$fits = main_est$fits
    if (undersmooth){
      results$fits_us = main_est$fits_us
    }
  }
  if (save_density_est){
    results$density_est = main_est$density_est
    results$density_preds = main_est$density_preds
  }
  return(results)
}
