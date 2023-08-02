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
#' @param mrk This indicates whether we want to estimate a multidimensional regression kink design.
#' @param kink This indicates whether we want to estimate a derivative discontinuity at the threshold. This option is typically used when we don't have an MRK/MRDK design but still want to estimate a discontinuity in the derivative.
#' @param eps This value indicates the step size used for the numerical derivative computation in the MRK estimation.
#' @param R1 The input for this option should be a vector containing indices for observations in treated region. Typically, this can be left blank, in which case R1 defaults to NULL.
#' @param R2 The input for this option should be a vector containing indices for observations in untreated region. Typically, this can be left blank, in which case R2 defaults to NULL.
#' @param K This is a parameter controlling the penalization for thin plate regression splines. Typically this is left empty, in which case K defaults to -1.
#' @param fuzzy The input should be a logical value (TRUE or FALSE) indicating whether this is a fuzzy RD. If this is not specified, the option defaults to FALSE (so the code assumes that this is a sharp RD).
#' @param reverse The input should be a logical value (TRUE or FALSE) indicating whether the treatment is reversed, i.e., what is defined as the "treatment group" for purposes of the MRD due to the region of the running variable space that defines it is is actually the control for the empirical application.
#' @param save_fits The input should be a logical value indicating whether the fitted thin plate regression splines should be saved. If the option is not specified, it defaults to TRUE.
#' @param save_density_est The input should be a logical value indicating whether kernel density estimates should be saved. If the option is not specified, it defaults to TRUE.
#' @param pred_dfs The input should be a dataframe with the values of the running variables that conditional average treatment effects should be calculated at. If left empty, it defaults to NULL.
#' @param bc The input should be a logical value indicating whether the thin plate regression splines should be undersmoothed by using the MSE-optimal penalty parameter from a higher-order spline as the penalty parameter for the splines actually used to compute the treatment effects.
#' @param undersmooth The input should be a logical value indicating whether thin plate regression splines should be undersmoothed by dividing the MSE-optimal penalty parameter by two, in order to reduce bias (for inference purposes).
#' @param save_sp The input should be a logical value indicating whether the parameter for thin plate regression spline fitting should be saved. If not specified, it defaults to TRUE.
#' @keywords basic_mrd_est
#' @export

basic_est = function(Y, X1, X2, W = NULL, F1.range = NULL, F2.range = NULL,
                     mrk, kink, eps,
                     R1 = NULL, R2 = NULL, K = -1,
                     fuzzy = FALSE, reverse = FALSE,
                     save_fits = TRUE, save_density_est = TRUE,
                     pred_dfs = NULL,
                     bc, undersmooth, save_sp = TRUE){
  ### Set F1 and F2 if they were not specified
  if (mrk){
    fuzzy = 1
  }
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
  fits = list(reduced = list(), stage1 = list())
  if (fuzzy){
    E = 3
  }else{
    E = 1
  }
  if (bc){
    fits_bc = fits
  }
  if (undersmooth){
    fits_us = fits
  }
  
  # Fit splines
  for (r in 1:2){
    fits$reduced[[r]] = gam(Y ~ s(X1, X2, k = K), subset = regions[[r]])
    if (bc){
      fit_higher_order = gam(Y ~ s(X1, X2, k = K, m=3), subset = regions[[r]])
      fits_bc$reduced[[r]] = gam(Y ~ s(X1, X2, k = K), subset = regions[[r]],
                                 sp = fit_higher_order$sp)
    }
    if (undersmooth){
      fits_us$reduced[[r]] = gam(Y ~ s(X1, X2, k = K), subset = regions[[r]],
                                 sp = 0.5*fits$reduced[[r]]$sp)
    }
    if (fuzzy){
      fits$stage1[[r]] = gam(W ~ s(X1, X2, k = K), subset = regions[[r]])
      if (bc){
        fit_higher_order = gam(W ~ s(X1, X2, k = K, m=3), subset = regions[[r]])
        fits_bc$stage1[[r]] = gam(W ~ s(X1, X2, k = K), subset = regions[[r]],
                                  sp = fit_higher_order$sp)
      }
      if (undersmooth){
        fits_us$stage1[[r]] = gam(W ~ s(X1, X2, k = K), subset = regions[[r]],
                                  sp = 0.5*fits$stage1[[r]]$sp)
      }
    }
  }
  
  # Obtain predictions
  preds = list(reduced = list(F1 = list(), F2 = list()),
               stage1 = list(F1 = list(), F2 = list()))
  preds_mrk = list(reduced = list(F1 = list(R1 = c(), R2 = c()), F2 = list(R1 = c(), R2 = c())),
                   stage1 = list(F1 = list(R1 = c(), R2 = c()), F2 = list(R1 = c(), R2 = c())))
  preds_bc = preds
  preds_mrk_bc = preds_mrk
  preds_us = preds
  preds_mrk_us = preds_mrk
  
  for (e in 1:max(E-1,1)){
    for (f in 1:2){
      for (r in 1:2){
        preds[[e]][[f]][[r]] = predict(fits[[e]][[r]], pred_dfs[[f]], se.fit = TRUE)
        if (bc){
          preds_bc[[e]][[f]][[r]] = predict(fits_bc[[e]][[r]], pred_dfs[[f]], se.fit = TRUE)
        }
        if (undersmooth){
          preds_us[[e]][[f]][[r]] = predict(fits_us[[e]][[r]], pred_dfs[[f]], se.fit = TRUE)
        }
        if (mrk | kink){
          # Create prediction data frames for numerical derivative estimation
          if (f == 1){
            pred_df_new1 = pred_dfs[[f]]
            pred_df_new0 = pred_dfs[[f]]
            if (r == 1){
              pred_df_new1$X1 = eps
              pred_df_new0$X1 = 0
            }else{
              pred_df_new1$X1 = 0
              pred_df_new0$X1 = -eps
            }
          }else{
            pred_df_new1 = pred_dfs[[f]]
            pred_df_new0 = pred_dfs[[f]]
            if (r == 1){
              pred_df_new1$X2 = eps
              pred_df_new0$X2 = 0
            }else{
              pred_df_new1$X2 = 0
              pred_df_new0$X2 = -eps
            }
          }
          # Compute numerical derivatives
          if (f == 1){
            pred1 = predict(fits[[e]][[r]],pred_df_new1,type="lpmatrix")
            pred0 = predict(fits[[e]][[r]],pred_df_new0,type="lpmatrix")
            preds_mrk[[e]][[f]][[r]]$basis_partial_X1 = (pred1-pred0)/eps
            preds_mrk[[e]][[f]][[r]]$partial_X1 = ((pred1-pred0)/eps) %*% fits[[e]][[r]]$coefficients
            if (bc){
              pred1_bc = predict(fits_bc[[e]][[r]],pred_df_new1,type="lpmatrix")
              pred0_bc = predict(fits_bc[[e]][[r]],pred_df_new0,type="lpmatrix")
              preds_mrk_bc[[e]][[f]][[r]]$basis_partial_X1 = (pred1_bc-pred0_bc)/eps
              preds_mrk_bc[[e]][[f]][[r]]$partial_X1 = ((pred1_bc-pred0_bc)/eps) %*% fits_bc[[e]][[r]]$coefficients
            }
            if (undersmooth){
              pred1_us = predict(fits_us[[e]][[r]],pred_df_new1,type="lpmatrix")
              pred0_us = predict(fits_us[[e]][[r]],pred_df_new0,type="lpmatrix")
              preds_mrk_us[[e]][[f]][[r]]$basis_partial_X1 = (pred1_us-pred0_us)/eps
              preds_mrk_us[[e]][[f]][[r]]$partial_X1 = ((pred1_us-pred0_us)/eps) %*% fits_us[[e]][[r]]$coefficients
            }
          }else{
            pred1 = predict(fits[[e]][[r]],pred_df_new1,type="lpmatrix")
            pred0 = predict(fits[[e]][[r]],pred_df_new0,type="lpmatrix")
            preds_mrk[[e]][[f]][[r]]$basis_partial_X2 = (pred1-pred0)/eps
            preds_mrk[[e]][[f]][[r]]$partial_X2 = ((pred1-pred0)/eps) %*% fits[[e]][[r]]$coefficients
            if (bc){
              pred1_bc = predict(fits_bc[[e]][[r]],pred_df_new1,type="lpmatrix")
              pred0_bc = predict(fits_bc[[e]][[r]],pred_df_new0,type="lpmatrix")
              preds_mrk_bc[[e]][[f]][[r]]$basis_partial_X2 = (pred1_bc-pred0_bc)/eps
              preds_mrk_bc[[e]][[f]][[r]]$partial_X2 = ((pred1_bc-pred0_bc)/eps) %*% fits_bc[[e]][[r]]$coefficients
            }
            if (undersmooth){
              pred1_us = predict(fits_us[[e]][[r]],pred_df_new1,type="lpmatrix")
              pred0_us = predict(fits_us[[e]][[r]],pred_df_new0,type="lpmatrix")
              preds_mrk_us[[e]][[f]][[r]]$basis_partial_X2 = (pred1_us-pred0_us)/eps
              preds_mrk_us[[e]][[f]][[r]]$partial_X2 = ((pred1_us-pred0_us)/eps) %*% fits_us[[e]][[r]]$coefficients
            }
          }
        }
      }
    }
  }
  # Treatment effects
  tau = list(reduced = list(pw = list(F1 = list(), F2 = list()),
                            agg = list(F1 = list(), F2 = list(), all = list())),
             stage1 = list(pw = list(F1 = list(), F2 = list()),
                           agg = list(F1 = list(), F2 = list(), all = list())),
             late = list(pw = list(F1 = list(), F2 = list()),
                         agg = list(F1 = list(), F2 = list(), all = list())))
  tau_bc = tau
  tau_us = tau
  tau_mrk = tau
  tau_mrk_bc = tau
  tau_mrk_us = tau
  
  # Compute pointwise treatment effects (CATE Estimates)
  for (f in 1:2){
    for (e in 1:max(E-1,1)){
      tau[[e]]$pw[[f]] = preds[[e]][[f]][[1]]$fit - preds[[e]][[f]][[2]]$fit
      if (bc){
        tau_bc[[e]]$pw[[f]] = preds_bc[[e]][[f]][[1]]$fit - preds_bc[[e]][[f]][[2]]$fit
      }
      if (undersmooth){
        tau_us[[e]]$pw[[f]] = preds_us[[e]][[f]][[1]]$fit - preds_us[[e]][[f]][[2]]$fit
      }
      if (mrk | kink){
        if (f == 1){
          tau_mrk[[e]]$pw[[f]] = preds_mrk[[e]][[f]][[1]]$partial_X1 - preds_mrk[[e]][[f]][[2]]$partial_X1
          if (bc){
            tau_mrk_bc[[e]]$pw[[f]] = preds_mrk_bc[[e]][[f]][[1]]$partial_X1 - preds_mrk_bc[[e]][[f]][[2]]$partial_X1
          }
          if (undersmooth){
            tau_mrk_us[[e]]$pw[[f]] = preds_mrk_us[[e]][[f]][[1]]$partial_X1 - preds_mrk_us[[e]][[f]][[2]]$partial_X1
          }
        }else{
          tau_mrk[[e]]$pw[[f]] = preds_mrk[[e]][[f]][[1]]$partial_X2 - preds_mrk[[e]][[f]][[2]]$partial_X2
          if (bc){
            tau_mrk_bc[[e]]$pw[[f]] = preds_mrk_bc[[e]][[f]][[1]]$partial_X2 - preds_mrk_bc[[e]][[f]][[2]]$partial_X2
          }
          if (undersmooth){
            tau_mrk_us[[e]]$pw[[f]] = preds_mrk_us[[e]][[f]][[1]]$partial_X2 - preds_mrk_us[[e]][[f]][[2]]$partial_X2
          }
        }
      }
    }
    # Reverse the estimates if what is defined as the "treatment group" for purposes of the MRD due to the region of the running variable space that defines it is is actually the control for the empirical application
    if (reverse){
      tau[[e]]$pw[[f]] = -tau[[e]]$pw[[f]]
      if (bc){
        tau_bc[[e]]$pw[[f]] = -tau_bc[[e]]$pw[[f]]
      }
      if (undersmooth){
        tau_us[[e]]$pw[[f]] = -tau_us[[e]]$pw[[f]]
      }
      if (mrk | kink){
        tau_mrk[[e]]$pw[[f]] = -tau_mrk[[e]]$pw[[f]]
        if (bc){
          tau_mrk_bc[[e]]$pw[[f]] = -tau_mrk_bc[[e]]$pw[[f]]
        }
        if (undersmooth){
          tau_mrk_us[[e]]$pw[[f]] = -tau_mrk_us[[e]]$pw[[f]]
        }
      }
    }
    # Compute the Wald ratio for the fuzzy design
    if (fuzzy){
      tau$late$pw[[f]] = tau$reduced$pw[[f]]/tau$stage1$pw[[f]]
      if (bc){
        tau_bc$late$pw[[f]] = tau_bc$reduced$pw[[f]]/tau_bc$stage1$pw[[f]]
      }
      if (undersmooth){
        tau_us$late$pw[[f]] = tau_us$reduced$pw[[f]]/tau_us$stage1$pw[[f]]
      }
    }
    if (mrk){
      tau_mrk$late$pw[[f]] = tau_mrk$reduced$pw[[f]]/tau_mrk$stage1$pw[[f]]
      if (bc){
        tau_mrk_bc$late$pw[[f]] = tau_mrk_bc$reduced$pw[[f]]/tau_mrk_bc$stage1$pw[[f]]
      }
      if (undersmooth){
        tau_mrk_us$late$pw[[f]] = tau_mrk_us$reduced$pw[[f]]/tau_mrk_us$stage1$pw[[f]]
      }
    }
    # Aggregate CATE estimates over the two segments of the treatment frontier
    for (e in 1:E){
      tau[[e]]$agg[[f]] = sum(tau[[e]]$pw[[f]]*kernel_est$mass[[f]])
      if (mrk | kink){
        tau_mrk[[e]]$agg[[f]] = sum(tau_mrk[[e]]$pw[[f]]*kernel_est$mass[[f]])
      }
      if (bc){
        tau_bc[[e]]$agg[[f]] = sum(tau_bc[[e]]$pw[[f]]*kernel_est$mass[[f]])
        if (mrk | kink){
          tau_mrk_bc[[e]]$agg[[f]] = sum(tau_mrk_bc[[e]]$pw[[f]]*kernel_est$mass[[f]])
        }
      }
      if (undersmooth){
        tau_us[[e]]$agg[[f]] = sum(tau_us[[e]]$pw[[f]]*kernel_est$mass[[f]])
        if (mrk | kink){
          tau_mrk_us[[e]]$agg[[f]] = sum(tau_mrk_us[[e]]$pw[[f]]*kernel_est$mass[[f]])
        }
      }
    }
  }
  
  # Aggregate CATE estimates over the entire treatment frontier
  for (e in 1:E){
    tau[[e]]$agg[[3]] = sum(c(tau[[e]]$pw[[1]],tau[[e]]$pw[[2]])*
                              kernel_est$mass[[3]])
    if (mrk | kink){
      tau_mrk[[e]]$agg[[3]] = sum(c(tau_mrk[[e]]$pw$F1,tau_mrk[[e]]$pw$F2)*kernel_est$mass[[3]])
    }
    if (bc){
      tau_bc[[e]]$agg[[3]] = sum(c(tau_bc[[e]]$pw[[1]],tau_bc[[e]]$pw[[2]])*
                                   kernel_est$mass[[3]])
      if (mrk | kink){
        tau_mrk_bc[[e]]$agg[[3]] = sum(c(tau_mrk_bc[[e]]$pw$F1,tau_mrk_bc[[e]]$pw$F2)*
                                         kernel_est$mass[[3]])
      }
    }
    if (undersmooth){
      tau_us[[e]]$agg[[3]] = sum(c(tau_us[[e]]$pw[[1]],tau_us[[e]]$pw[[2]])*
                                   kernel_est$mass[[3]])
      if (mrk | kink){
        tau_mrk_us[[e]]$agg[[3]] = sum(c(tau_mrk_us[[e]]$pw$F1,tau_mrk_us[[e]]$pw$F2)*
                                         kernel_est$mass[[3]])
      }
    }
  }
  # Save Bayesian standard errors if desired
  tau_se = list(reduced = list(pw = list(bayes = list(list(F1=list(), F2=list()))),
                               agg=list(delta = list(F1=list(), F2=list(), all=list()))),
                stage1 = list(pw = list(bayes = list(list(F1=list(), F2=list()))),
                              agg=list(delta = list(F1=list(), F2=list(), all=list()))))
  tau_se_bc = tau_se
  tau_se_us = tau_se
  
  for (e in 1:max(E-1,1)){
    for (f in 1:2){
      tau_se[[e]]$pw$bayes[[f]] =
        sqrt(preds[[e]][[f]][[1]]$se.fit^2 + preds[[e]][[f]][[2]]$se.fit^2)
      if (bc){
        tau_se_bc[[e]]$pw$bayes[[f]] =
          sqrt(preds_bc[[e]][[f]][[1]]$se.fit^2 + preds_bc[[e]][[f]][[2]]$se.fit^2)
      }
      if (undersmooth){
        tau_se_us[[e]]$pw$bayes[[f]] =
          sqrt(preds_us[[e]][[f]][[1]]$se.fit^2 + preds_us[[e]][[f]][[2]]$se.fit^2)
      }
    }
  }
  results = list(tau = tau, tau_bc = tau_bc, tau_us = tau_us,
                 tau_mrk = tau_mrk, tau_mrk_bc = tau_mrk_bc, tau_mrk_us = tau_mrk_us,
                 preds = preds, preds_bc = preds_bc, preds_us = preds_us,
                 preds_mrk = preds_mrk, preds_mrk_bc = preds_mrk_bc, preds_mrk_us = preds_mrk_us,
                 tau_se = tau_se, tau_se_bc = tau_se_bc, tau_se_us = tau_se_us)
  if (save_fits){
    results$fits = fits
    if (bc){
      results$fits_bc = fits_bc
    }
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
    if (bc){
      results$sp_bc = list(reduced = c(fits_bc$reduced[[1]]$sp,fits_bc$reduced[[2]]$sp))
    }
    if (undersmooth){
      results$sp_us = list(reduced = c(fits_us$reduced[[1]]$sp,fits_us$reduced[[2]]$sp))
    }
    if (fuzzy){
      results$sp$stage1 = c(fits$stage1[[1]]$sp,fits$stage1[[2]]$sp)
      if (bc){
        results$sp_bc$stage1 = c(fits_bc$stage1[[1]]$sp,fits_bc$stage1[[2]]$sp)
      }
      if (undersmooth){
        results$sp_us$stage1 = c(fits_us$stage1[[1]]$sp,fits_us$stage1[[2]]$sp)
      }
    }
  }
  return(results)
}

#' Main function for MRD estimation.
#' 
#' This function carries out estimation for regression discontinuity designs with two running variables (MRD) and regression kink designs with two running variables (MRK).
#' The function is written assuming that the treatment frontier (i.e., the boundary in running variable space dividing the treated and untreated regions) is given by the positive X1 and X2 axes. In some empirical applications where these are not the exact conditions for program eligibility (e.g., the threshold is not zero, a running variable must fall below instead of above a threshold, or observations in the positive quadrant are untreated), a simple redefinition of the running variables (e.g., through translation or by multiplying by negative one) or a redefinition of treatment (to label treated units as untreated, and then multiplying the final estimates by negative one) is all that is necessary. In cases where the treatment frontier is still a boundary separating treated and untreated regions but is not given by two straight lines forming a right angle, this function can still be used in principle (by explicitly specifying the indices for treated and untreated regions via "R1" and "R2", and by specifying the coordinates of the treatment frontier via "pred_dfs"), but this functionality has not been extensively tested.
#' Please feel free to contact the author with questions at alden15@nber.org.
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
#' @param mrk This indicates whether we want to estimate a multidimensional regression kink design.
#' @param kink This indicates whether we want to estimate a derivative discontinuity at the threshold. This option is typically used when we don't have an MRK/MRDK design but still want to estimate a discontinuity in the derivative.
#' @param sharp_mrk This indicates whether this is a sharp MRK design. If this is set to TRUE, then analytic standard errors for the MRK estimate are computed ignoring the estimation uncertainty in the denominator.
#' @param sharp_mrd This indicates whether this is a sharp MRD design. If this is set to TRUE but the denominator is also estimated, then analytic standard errors for the MRD estimate are computed ignoring the estimation uncertainty in the denominator.
#' @param eps This value indicates the step size used for the numerical derivative computation in the MRK estimation.
#' @param R1 The input for this option should be a vector containing indices for observations in treated region. Typically, this can be left blank, in which case R1 defaults to NULL.
#' @param R2 The input for this option should be a vector containing indices for observations in untreated region. Typically, this can be left blank, in which case R2 defaults to NULL.
#' @param signif_lvl The indicates the statistical significance level for computing the simultaneous confidence band. If left empty, this defaults to the 5 percent significance level.
#' @param K This is a parameter controlling the penalization for thin plate regression splines. Typically this is left empty, in which case K defaults to -1.
#' @param bs The input should be a logical value indicating whether to use bootstrap for standard errors. If not specified, this option defaults to FALSE.
#' @param bs_reps The input should be the number of desired bootstrap iterations. If not specified, this is option defaults to 1000.
#' @param bs_seeds The input should be a logical value indicating whether the random seed for the rth bootstrap iteration should be set to r. If not specified, this option defaults to TRUE.
#' @param save_density_est The input should be a logical value indicating whether kernel density estimates should be saved. If the option is not specified, it defaults to TRUE.
#' @param save_fits The input should be a logical value indicating whether the fitted thin plate regression splines should be saved. If the option is not specified, it defaults to TRUE.
#' @param fuzzy The input should be a logical value (TRUE or FALSE) indicating whether this is a fuzzy RD. If this is not specified, the option defaults to FALSE (so the code assumes that this is a sharp RD).
#' @param reverse The input should be a logical value (TRUE or FALSE) indicating whether the treatment is reversed, i.e., what is defined as the "treatment group" for purposes of the MRD due to the region of the running variable space that defines it is is actually the control for the empirical application.
#' @param fn_ci This indicates whether simultaneous confidence bands should be computed for the CATE estimates along the treatment frontier.
#' @param reps_fn_ci This indicates the number of simulation draws to be used for computing the simultaneous confidence bands for the CATE estimates.
#' @param bc The input should be a logical value indicating whether the thin plate regression splines should be undersmoothed by using the MSE-optimal penalty parameter from a higher-order spline as the penalty parameter for the splines actually used to compute the treatment effects.
#' @param undersmooth The input should be a logical value indicating whether thin plate regression splines should be undersmoothed by dividing the MSE-optimal penalty parameter by two, in order to reduce bias (for inference purposes).
#' @param pred_dfs The input should be a dataframe with the values of the running variables that conditional average treatment effects should be calculated at. If left empty, it defaults to NULL.
#' @param bs_print This indicates the multiples at which bootstrap iterations should be printed (if the bs option is set to TRUE).
#' @return The mrdd.fn returns a nested list containing the estimation results. Key components of the returned list include:
#' \itemize{
#'  \item{"tau"}{This sublist includes the MRD estimates. Furthered nested within this sublist may include the subsublists: "reduced", "stage1", and "late". For sharp MRD designs, only "reduced" is populated and contains the MRD estimate, whereas for fuzzy designs, "stage1" contains estimates for the denominator of the Wald ratio, and "late" is the local average treatment effect/Wald Ratio. Further nested within "reduced" (as well as "stage1" and "late") are the subsubsublists "pw" and "agg", which contain either the pointwise CATE estimates, or the CATE estimates aggregated over certain parts of the treatment frontier respectively.}
#'  \item{"tau_se"}{This sublist contains standard errors for the MRD estimates. It contains similar subsublists to the "tau" sublist, which in turn contains standard errors for the CATE at every point along the treatment frontier ("pw"), or for the aggregated treatment effect estimate ("agg"). Nested within "pw" and "delta" may include "bayes" which corresponds to the Bayesian standard errors (for "pw") or "delta" which corresponds to the delta method (for "agg") which incorporates both uncertainty in the thin plate regression spline estimation and the density estimation, and if the "bs" option in the function is specified to be true, then bootstrap standard errors ("bs") are also available.}
#'  \item{"fn_cov"}{This sublist contains critical values for the simultaneous confidence bands (corresponding to the significance level as specified by the "signif_lvl" option in the mrdd.fn function). For example, while the critical value for constructing a 95 percent pointwise confidence interval is 1.96, it is typically larger for a 95 percent simultaneous confidence band, which should contain the entire CATE over the treatment frontier 95 percent of the time over random samples.}
#'  \item{"fits"}{This sublists contains the fitted thin plate regression splines underlying the MRD/MRK estimates.}
#'. \item{"density_est" and "density_preds"}{These sublists contain density estimation results and the predicted values of the density function over the treatment frontier, which are used for computing the aggregated treatment effect over the treatment frontier.}
#' }
#' Many variations of the sublists mentioned above are included in the brief description above are also included in the output. To better understand what some of them mean, the suffix "bc" and "us" refer to quantities associated with estimation that uses bias-corrected estimates (based on higher-order thin plate regression splines), or undersmoothed estimates (by dividing the MSE-optimal penalty parameter by one-half) respectively. The suffix "mrk" refers to multidimensional regression kink estimates (or if the "mrk" option is FALSE and the "kink" option is TRUE, simply a kink estimate), The terms "F1" and "F2" refer to the region of the treatment frontier where the first running variable is equal to zero and the second is positive, and vice versa, respectively.
#' The terms "reduced", "stage1", and "late" are borrowed from the instrumental variables literature. For sharp MRD designs, the MRD estimate corresponds to the results under "reduced". For fuzzy designs (and MRK designs), "reduced" and "stage1" correspond to the numerator and denominator of the Wald ratio, and the MRD (or MRK) treatment effect is given by the Wald ratio itself, abbreviated as "late" (short for local average treatment effect).

mrdd.fn = function(Y, X1, X2, W = NULL, F1.range = NULL, F2.range = NULL,
                   mrk = FALSE, kink = TRUE, sharp_mrk = TRUE, sharp_mrd = FALSE, eps = 10^(-7),
                   R1 = NULL, R2 = NULL,
                   signif_lvl = 0.05, K = -1,
                   bs = FALSE, bs_reps = 1000, bs_seeds = TRUE,
                   save_density_est = TRUE, save_fits = TRUE,
                   fuzzy = FALSE, reverse = FALSE, fn_ci = TRUE, reps_fn_ci = 1000,
                   bc = TRUE, undersmooth = TRUE,
                   pred_dfs = NULL, bs_print = 100){
  if (mrk){
    fuzzy = 1
  }
  if (!mrk){
    sharp_mrk = FALSE
  }
  zscore = qnorm(1-signif_lvl/2)
  if (is.null(F1.range)){
    F1.range = seq(0, max(X2), length = 100)
  }
  if (is.null(F2.range)){
    F2.range = seq(0, max(X1), length = 100)
  }
  # Main estimates
  main_est = basic_est(Y = Y, X1 = X1, X2 = X2, W = W, K = K,
                       mrk = mrk, kink = kink, eps=eps,
                       F1.range = F1.range, F2.range = F2.range, R1 = R1, R2 = R2,
                       save_density_est = save_density_est,
                       fuzzy=fuzzy, reverse = reverse,
                       bc = bc, undersmooth = undersmooth, save_sp = TRUE,
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
  
  # Create lists/dataframes to store bootstrap results
  if (is.null(pred_dfs)){
    bs_pred_dfs = list(data.frame(bs_X1=0, bs_X2=F1.range),
                       data.frame(bs_X2=0, bs_X1=F2.range))
  }else{
    bs_pred_dfs = list(data.frame(bs_X1 = pred_dfs[[1]]$X1, bs_X2 = pred_dfs[[1]]$X2),
                       data.frame(bs_X1 = pred_dfs[[2]]$X1, bs_X2 = pred_dfs[[2]]$X2))
  }
  bs_est = list(reduced = list(pw = list(F1 = list(), F2 = list()),
                               agg = list(F1 = list(), F2 = list(), all = list())),
                stage1 = list(pw = list(F1 = list(), F2 = list()),
                              agg = list(F1 = list(), F2 = list(), all = list())),
                late = list(pw = list(F1 = list(), F2 = list()),
                            agg = list(F1 = list(), F2 = list(), all = list())))
  bs_est$density_est = list(F1 = matrix(nrow = length(F1.range), ncol = bs_reps),
                            F2 = matrix(nrow = length(F2.range), ncol = bs_reps),
                            all = matrix(nrow = length(F1.range) + length(F2.range), ncol = bs_reps))
  bs_est$density_preds = list(F1 = matrix(nrow = length(F1.range), ncol = bs_reps),
                              F2 = matrix(nrow = length(F2.range), ncol = bs_reps),
                              all = matrix(nrow = length(F1.range) + length(F2.range), ncol = bs_reps))
  bs_est_bc = bs_est
  bs_est_us = bs_est
  bs_est_mrk = bs_est
  bs_est_mrk_bc = bs_est
  bs_est_mrk_us = bs_est
  
  for (e in 1:E){
    for (f in 1:2){
      bs_est[[e]]$pw[[f]] =
        matrix(nrow = length(main_est$tau$reduced$pw[[f]]), ncol = bs_reps)
      bs_est[[e]]$agg[[f]] = rep(NA, bs_reps)
      if (bc){
        bs_est_bc[[e]]$pw[[f]] =
          matrix(nrow = length(main_est$tau_bc$reduced$pw[[f]]), ncol = bs_reps)
        bs_est_bc[[e]]$agg[[f]] = rep(NA, bs_reps)
      }
      if (undersmooth){
        bs_est_us[[e]]$pw[[f]] =
          matrix(nrow = length(main_est$tau_us$reduced$pw[[f]]), ncol = bs_reps)
        bs_est_us[[e]]$agg[[f]] = rep(NA, bs_reps)
      }
      if (mrk | kink){
        bs_est_mrk[[e]]$pw[[f]] =
          matrix(nrow = length(main_est$tau_mrk$reduced$pw[[f]]), ncol = bs_reps)
        bs_est_mrk[[e]]$agg[[f]] = rep(NA, bs_reps)
        if (bc){
          bs_est_mrk_bc[[e]]$pw[[f]] =
            matrix(nrow = length(main_est$tau_mrk_bc$reduced$pw[[f]]), ncol = bs_reps)
          bs_est_mrk_bc[[e]]$agg[[f]] = rep(NA, bs_reps)
        }
        if (undersmooth){
          bs_est_mrk_us[[e]]$pw[[f]] =
            matrix(nrow = length(main_est$tau_mrk_us$reduced$pw[[f]]), ncol = bs_reps)
          bs_est_mrk_us[[e]]$agg[[f]] = rep(NA, bs_reps)
        }
      }
    }
    bs_est[[e]]$agg[[3]] = rep(NA, bs_reps)
    if (bc){
      bs_est_bc[[e]]$agg[[3]] = rep(NA, bs_reps)
    }
    if (undersmooth){
      bs_est_us[[e]]$agg[[3]] = rep(NA, bs_reps)
    }
    if (mrk | kink){
      bs_est_mrk[[e]]$agg[[3]] = rep(NA, bs_reps)
      if (bc){
        bs_est_mrk_bc[[e]]$agg[[3]] = rep(NA, bs_reps)
      }
      if (undersmooth){
        bs_est_mrk_us[[e]]$agg[[3]] = rep(NA, bs_reps)
      }
    }
  }
  
  # Do bootstrap (only for kernel density estimation if bs option not selected)
  for (b in 1:bs_reps){
    if (bs_seeds){
      set.seed(b)
    }
    bs_indices = sample.int(n = length(Y), size = length(Y), replace = TRUE)
    bs_X1 = X1[bs_indices]
    bs_X2 = X2[bs_indices]
    
    # Bootstrap for KDE
    kernel_est = list()
    kernel_density_est = locfit(~ bs_X1 + bs_X2)
    for (f in 1:2){
      kernel_est$preds[[f]] = predict(kernel_density_est, bs_pred_dfs[[f]])
      kernel_est$mass[[f]] = kernel_est$preds[[f]]/sum(kernel_est$preds[[f]])
    }
    kernel_est$mass[[3]] = c(kernel_est$preds[[1]],kernel_est$preds[[2]])/
      sum(c(kernel_est$preds[[1]],kernel_est$preds[[2]]))
    
    for (f in 1:3){
      bs_est$density_est[[f]][,b] = kernel_est$mass[[f]]
      bs_est_us$density_est[[f]][,b] = kernel_est$mass[[f]]
      bs_est_mrk$density_est[[f]][,b] = kernel_est$mass[[f]]
      bs_est_mrk_us$density_est[[f]][,b] = kernel_est$mass[[f]]
      if (f != 3){
        bs_est$density_preds[[f]][,b] = kernel_est$preds[[f]]
        bs_est_us$density_preds[[f]][,b] = kernel_est$preds[[f]]
        bs_est_mrk$density_preds[[f]][,b] = kernel_est$preds[[f]]
        bs_est_mrk_us$density_preds[[f]][,b] = kernel_est$preds[[f]]
      }
    }
    # Bootstrap for MRD/MRK estimation
    if (bs){
      if (!(b%%bs_print)){
        print(paste0("Bootstrap Iteration ", b, " of ", bs_reps))
      }
      bs_Y = Y[bs_indices]
      if (!is.null(W)){
        bs_W = W[bs_indices]
      }
      current_bs_est = basic_est(Y = bs_Y, X1 = bs_X1, X2 = bs_X2, W = bs_W, K=K,
                                 mrk = mrk, kink = kink, eps=eps,
                                 F1.range = F1.range, F2.range = F2.range, R1 = R1, R2 = R2,
                                 save_density_est = save_density_est,
                                 fuzzy=fuzzy, reverse = reverse,
                                 bc = bc, undersmooth = undersmooth, save_sp = TRUE,
                                 pred_dfs = pred_dfs)
      for (e in 1:E){
        for (f in 1:2){
          bs_est[[e]]$pw[[f]][,b] = current_bs_est$tau[[e]]$pw[[f]]
          bs_est[[e]]$agg[[f]][b] = current_bs_est$tau[[e]]$agg[[f]]
          if (bc){
            bs_est_bc[[e]]$pw[[f]][,b] = current_bs_est$tau_bc[[e]]$pw[[f]]
            bs_est_bc[[e]]$agg[[f]][b] = current_bs_est$tau_bc[[e]]$agg[[f]]
          }
          if (undersmooth){
            bs_est_us[[e]]$pw[[f]][,b] = current_bs_est$tau_us[[e]]$pw[[f]]
            bs_est_us[[e]]$agg[[f]][b] = current_bs_est$tau_us[[e]]$agg[[f]]
          }
          if (mrk | kink){
            bs_est_mrk[[e]]$pw[[f]][,b] = current_bs_est$tau_mrk[[e]]$pw[[f]]
            bs_est_mrk[[e]]$agg[[f]][b] = current_bs_est$tau_mrk[[e]]$agg[[f]]
            if (bc){
              bs_est_mrk_bc[[e]]$pw[[f]][,b] = current_bs_est$tau_mrk_bc[[e]]$pw[[f]]
              bs_est_mrk_bc[[e]]$agg[[f]][b] = current_bs_est$tau_mrk_bc[[e]]$agg[[f]]
            }
            if (undersmooth){
              bs_est_mrk_us[[e]]$pw[[f]][,b] = current_bs_est$tau_mrk_us[[e]]$pw[[f]]
              bs_est_mrk_us[[e]]$agg[[f]][b] = current_bs_est$tau_mrk_us[[e]]$agg[[f]]
            }
          }
        }
        bs_est[[e]]$agg[[3]][b] = current_bs_est$tau[[e]]$agg[[3]]
        if (bc){
          bs_est_bc[[e]]$agg[[3]][b] = current_bs_est$tau_bc[[e]]$agg[[3]]
        }
        if (undersmooth){
          bs_est_us[[e]]$agg[[3]][b] = current_bs_est$tau_us[[e]]$agg[[3]]
        }
        if (mrk | kink){
          bs_est_mrk[[e]]$agg[[3]][b] = current_bs_est$tau_mrk[[e]]$agg[[3]]
          if (bc){
            bs_est_mrk_bc[[e]]$agg[[3]][b] = current_bs_est$tau_mrk_bc[[e]]$agg[[3]]
          }
          if (undersmooth){
            bs_est_mrk_us[[e]]$agg[[3]][b] = current_bs_est$tau_mrk_us[[e]]$agg[[3]]
          }
        }
      }
    }
  }
  
  tau_se = main_est$tau_se
  tau_se_bc = main_est$tau_se_bc
  tau_se_us = main_est$tau_se_us
  
  if (fuzzy){
    tau_se$late = list(pw = list(), agg = list())
    if (bc){
      tau_se_bc$late = list(pw = list(), agg = list())
    }
    if (undersmooth){
      tau_se_us$late = list(pw = list(), agg = list())
    }
  }
  
  # Create objects to store variance estimates
  var_mat = list(reduced = list(pw = list(F1 = list(), F2 = list()), all = list()),
                 density_est = list(pw = list(F1 = list(), F2 = list(), all = list())))
  var_mat_bc = var_mat
  var_mat_us = var_mat
  var_mat_mrk = var_mat
  var_mat_mrk_bc = var_mat
  var_mat_mrk_us = var_mat
  
  tau_se_mrk = list(reduced = list(pw = list(F1 = list(), F2 = list()),
                                   agg = list(delta = list(F1 = list(), F2 = list(), all = list()))),
                    stage1 = list(pw = list(F1 = list(), F2 = list()),
                                  agg = list(delta = list(F1 = list(), F2 = list(), all = list()))),
                    late = list(pw = list(F1 = list(), F2 = list()),
                                agg = list(delta = list(F1 = list(), F2 = list(), all = list()))))
  tau_se_mrk_bc = tau_se_mrk
  tau_se_mrk_us = tau_se_mrk
  
  max_vec = list(reduced = list(pw = list(F1 = list(), F2 = list())),
                 stage1 = list(pw = list(F1 = list(), F2 = list())))
  max_vec_bc = max_vec
  max_vec_us = max_vec
  fn_ci_crit = list(reduced = list(F1 = list(), F2 = list()),
                    stage1 = list(F1 = list(), F2 = list()))
  fn_ci_crit_bc = fn_ci_crit
  fn_ci_crit_us = fn_ci_crit
  
  if (fuzzy){
    var_mat$stage1 = list(pw = list(F1 = list(), F2 = list()))
    var_mat_bc$stage1 = var_mat$stage1
    var_mat_us$stage1 = var_mat$stage1
    var_mat_mrk$stage1 = var_mat$stage1
    var_mat_mrk_bc$stage1 = var_mat$stage1
    var_mat_mrk_us$stage1 = var_mat$stage1
  }
  
  ### Compute standard errors for aggregated treatment effect estimates
  # Standard error for MSE-optimal estimates
  for (e in 1:max(E-1,1)){
    for (f in 1:2){
      # Compute covariance matrix of treatment effect point estimates
      var_mat[[e]][[f]] = 
        predict(main_est$fits[[e]][[1]], pred_dfs[[f]], type = "lpmatrix") %*% 
        main_est$fits[[e]][[1]]$Vp %*%
        t(predict(main_est$fits[[e]][[1]], pred_dfs[[f]], type = "lpmatrix")) + 
        predict(main_est$fits[[e]][[2]], pred_dfs[[f]], type = "lpmatrix") %*% 
        main_est$fits[[e]][[2]]$Vp %*%
        t(predict(main_est$fits[[e]][[2]], pred_dfs[[f]], type = "lpmatrix"))
      if (mrk | kink){
        if (f == 1){
          var_mat_mrk[[e]][[f]] = 
            main_est$preds_mrk[[e]][[f]][[1]]$basis_partial_X1 %*% 
            main_est$fits[[e]][[1]]$Vp %*%
            t(main_est$preds_mrk[[e]][[f]][[1]]$basis_partial_X1) + 
            main_est$preds_mrk[[e]][[f]][[2]]$basis_partial_X1 %*% 
            main_est$fits[[e]][[2]]$Vp %*%
            t(main_est$preds_mrk[[e]][[f]][[2]]$basis_partial_X1)
          tau_se_mrk[[e]]$pw[[f]] = sqrt(diag(var_mat_mrk[[e]][[f]]))
        }else{
          var_mat_mrk[[e]][[f]] = 
            main_est$preds_mrk[[e]][[f]][[1]]$basis_partial_X2 %*% 
            main_est$fits[[e]][[1]]$Vp %*%
            t(main_est$preds_mrk[[e]][[f]][[1]]$basis_partial_X2) + 
            main_est$preds_mrk[[e]][[f]][[2]]$basis_partial_X2 %*% 
            main_est$fits[[e]][[2]]$Vp %*%
            t(main_est$preds_mrk[[e]][[f]][[2]]$basis_partial_X2)
          tau_se_mrk[[e]]$pw[[f]] = sqrt(diag(var_mat_mrk[[e]][[f]]))
        }
      }
      # Compute covariance matrix for density estimate
      var_mat$density_est[[f]] = var(t(bs_est$density_preds[[f]]))/length(Y)
      # Compute standard error for aggregated treatment effect estimate using the delta method
      Deriv_wrt_f = as.numeric( (1/(sum(main_est$density_preds[[f]]))^2) * 
                                  ( sum(main_est$density_preds[[f]]) * main_est$tau[[e]]$pw[[f]] - 
                                      sum(main_est$density_preds[[f]] * main_est$tau[[e]]$pw[[f]]) ) )
      Deriv_wrt_tau = main_est$density_est[[f]]
      tau_se[[e]]$agg$delta[[f]] = sqrt(t(Deriv_wrt_tau) %*% var_mat[[e]][[f]] %*% Deriv_wrt_tau + 
                                          t(Deriv_wrt_f) %*% var_mat$density_est[[f]] %*% Deriv_wrt_f)
      # Compute standard error for aggregated MRK treatment effect estimate using the delta method
      if (mrk | kink){
        var_mat_mrk$density_est[[f]] = var_mat$density_est[[f]]
        Deriv_wrt_f_mrk = as.numeric( (1/(sum(main_est$density_preds[[f]]))^2) * 
                                        ( sum(main_est$density_preds[[f]]) * main_est$tau_mrk[[e]]$pw[[f]] - 
                                            sum(main_est$density_preds[[f]] * main_est$tau_mrk[[e]]$pw[[f]]) ) )
        Deriv_wrt_tau_mrk = main_est$density_est[[f]]
        tau_se_mrk[[e]]$agg$delta[[f]] = sqrt(t(Deriv_wrt_tau_mrk) %*% var_mat_mrk[[e]][[f]] %*% Deriv_wrt_tau_mrk + 
                                                t(Deriv_wrt_f_mrk) %*% var_mat_mrk$density_est[[f]] %*% Deriv_wrt_f_mrk)
      }
      
      # Standard error for undersmoothed estimates (using penalty parameter from higher-order spline)
      if (bc){
        # Compute covariance matrix of treatment effect point estimates
        var_mat_bc[[e]][[f]] = 
          predict(main_est$fits_bc[[e]][[1]], pred_dfs[[f]], type = "lpmatrix") %*% 
          main_est$fits_bc[[e]][[1]]$Vp %*%
          t(predict(main_est$fits_bc[[e]][[1]], pred_dfs[[f]], type = "lpmatrix")) + 
          predict(main_est$fits_bc[[e]][[2]], pred_dfs[[f]], type = "lpmatrix") %*% 
          main_est$fits_bc[[e]][[2]]$Vp %*%
          t(predict(main_est$fits_bc[[e]][[2]], pred_dfs[[f]], type = "lpmatrix"))
        if (mrk | kink){
          if (f == 1){
            var_mat_mrk_bc[[e]][[f]] = 
              main_est$preds_mrk_bc[[e]][[f]][[1]]$basis_partial_X1 %*% 
              main_est$fits_bc[[e]][[1]]$Vp %*%
              t(main_est$preds_mrk_bc[[e]][[f]][[1]]$basis_partial_X1) + 
              main_est$preds_mrk_bc[[e]][[f]][[2]]$basis_partial_X1 %*% 
              main_est$fits_bc[[e]][[2]]$Vp %*%
              t(main_est$preds_mrk_bc[[e]][[f]][[2]]$basis_partial_X1)
            tau_se_mrk_bc[[e]]$pw[[f]] = sqrt(diag(var_mat_mrk_bc[[e]][[f]]))
          }else{
            var_mat_mrk_bc[[e]][[f]] = 
              main_est$preds_mrk_bc[[e]][[f]][[1]]$basis_partial_X2 %*% 
              main_est$fits_bc[[e]][[1]]$Vp %*%
              t(main_est$preds_mrk_bc[[e]][[f]][[1]]$basis_partial_X2) + 
              main_est$preds_mrk_bc[[e]][[f]][[2]]$basis_partial_X2 %*% 
              main_est$fits_bc[[e]][[2]]$Vp %*%
              t(main_est$preds_mrk_bc[[e]][[f]][[2]]$basis_partial_X2)
            tau_se_mrk_bc[[e]]$pw[[f]] = sqrt(diag(var_mat_mrk_bc[[e]][[f]]))
          }
        }
        # Compute covariance matrix for density estimate
        var_mat_bc$density_est[[f]] = var_mat$density_est[[f]]
        # Compute standard error for aggregated treatment effect estimate using delta method
        Deriv_wrt_f = as.numeric( (1/(sum(main_est$density_preds[[f]]))^2) * 
                                    ( sum(main_est$density_preds[[f]]) * main_est$tau_bc[[e]]$pw[[f]] - 
                                        sum(main_est$density_preds[[f]] * main_est$tau_bc[[e]]$pw[[f]]) ) )
        Deriv_wrt_tau = main_est$density_est[[f]]
        tau_se_bc[[e]]$agg$delta[[f]] = sqrt(t(Deriv_wrt_tau) %*% var_mat_bc[[e]][[f]] %*% Deriv_wrt_tau + 
                                               t(Deriv_wrt_f) %*% var_mat_bc$density_est[[f]] %*% Deriv_wrt_f)
        # Compute standard error for aggregated MRK treatment effect estimate using delta method
        if (mrk | kink){
          var_mat_mrk_bc$density_est[[f]] = var_mat$density_est[[f]]
          Deriv_wrt_f_mrk_bc = as.numeric( (1/(sum(main_est$density_preds[[f]]))^2) * 
                                             ( sum(main_est$density_preds[[f]]) * main_est$tau_mrk_bc[[e]]$pw[[f]] - 
                                                 sum(main_est$density_preds[[f]] * main_est$tau_mrk_bc[[e]]$pw[[f]]) ) )
          Deriv_wrt_tau_mrk_bc = main_est$density_est[[f]]
          tau_se_mrk_bc[[e]]$agg$delta[[f]] = sqrt(t(Deriv_wrt_tau_mrk_bc) %*% var_mat_mrk_bc[[e]][[f]] %*% Deriv_wrt_tau_mrk_bc + 
                                                     t(Deriv_wrt_f_mrk_bc) %*% var_mat_mrk_bc$density_est[[f]] %*% Deriv_wrt_f_mrk_bc)
        }
      }
      
      # Standard error for undersmoothed estimates (using half of MSE-optimal penalty parameter)
      if (undersmooth){
        # Compute covariance matrix of treatment effect point estimates
        var_mat_us[[e]][[f]] = 
          predict(main_est$fits_us[[e]][[1]], pred_dfs[[f]], type = "lpmatrix") %*% 
          main_est$fits_us[[e]][[1]]$Vp %*%
          t(predict(main_est$fits_us[[e]][[1]], pred_dfs[[f]], type = "lpmatrix")) + 
          predict(main_est$fits_us[[e]][[2]], pred_dfs[[f]], type = "lpmatrix") %*% 
          main_est$fits_us[[e]][[2]]$Vp %*%
          t(predict(main_est$fits_us[[e]][[2]], pred_dfs[[f]], type = "lpmatrix"))
        if (mrk | kink){
          if (f == 1){
            var_mat_mrk_us[[e]][[f]] = 
              main_est$preds_mrk_us[[e]][[f]][[1]]$basis_partial_X1 %*% 
              main_est$fits_us[[e]][[1]]$Vp %*%
              t(main_est$preds_mrk_us[[e]][[f]][[1]]$basis_partial_X1) + 
              main_est$preds_mrk_us[[e]][[f]][[2]]$basis_partial_X1 %*% 
              main_est$fits_us[[e]][[2]]$Vp %*%
              t(main_est$preds_mrk_us[[e]][[f]][[2]]$basis_partial_X1)
            tau_se_mrk_us[[e]]$pw[[f]] = sqrt(diag(var_mat_mrk_us[[e]][[f]]))
          }else{
            var_mat_mrk_us[[e]][[f]] = 
              main_est$preds_mrk_us[[e]][[f]][[1]]$basis_partial_X2 %*% 
              main_est$fits_us[[e]][[1]]$Vp %*%
              t(main_est$preds_mrk_us[[e]][[f]][[1]]$basis_partial_X2) + 
              main_est$preds_mrk_us[[e]][[f]][[2]]$basis_partial_X2 %*% 
              main_est$fits_us[[e]][[2]]$Vp %*%
              t(main_est$preds_mrk_us[[e]][[f]][[2]]$basis_partial_X2)
            tau_se_mrk_us[[e]]$pw[[f]] = sqrt(diag(var_mat_mrk_us[[e]][[f]]))
          }
        }
        # Compute covariance matrix for density estimate
        var_mat_us$density_est[[f]] = var_mat$density_est[[f]]
        # Compute standard error for aggregated treatment effect estimate using delta method
        Deriv_wrt_f = as.numeric( (1/(sum(main_est$density_preds[[f]]))^2) * 
                                    ( sum(main_est$density_preds[[f]]) * main_est$tau_us[[e]]$pw[[f]] - 
                                        sum(main_est$density_preds[[f]] * main_est$tau_us[[e]]$pw[[f]]) ) )
        Deriv_wrt_tau = main_est$density_est[[f]]
        tau_se_us[[e]]$agg$delta[[f]] = sqrt(t(Deriv_wrt_tau) %*% var_mat_us[[e]][[f]] %*% Deriv_wrt_tau + 
                                               t(Deriv_wrt_f) %*% var_mat_us$density_est[[f]] %*% Deriv_wrt_f)
        # Compute standard error for aggregated MRK treatment effect estimate using delta method
        if (mrk | kink){
          var_mat_mrk_us$density_est[[f]] = var_mat$density_est[[f]]
          Deriv_wrt_f_mrk_us = as.numeric( (1/(sum(main_est$density_preds[[f]]))^2) * 
                                             ( sum(main_est$density_preds[[f]]) * main_est$tau_mrk_us[[e]]$pw[[f]] - 
                                                 sum(main_est$density_preds[[f]] * main_est$tau_mrk_us[[e]]$pw[[f]]) ) )
          Deriv_wrt_tau_mrk_us = main_est$density_est[[f]]
          tau_se_mrk_us[[e]]$agg$delta[[f]] = sqrt(t(Deriv_wrt_tau_mrk_us) %*% var_mat_mrk_us[[e]][[f]] %*% Deriv_wrt_tau_mrk_us + 
                                                     t(Deriv_wrt_f_mrk_us) %*% var_mat_mrk_us$density_est[[f]] %*% Deriv_wrt_f_mrk_us)
        }
      }
    }
    
    ### Compute standard error for CATE aggregated over the entire treatment frontier
    var_mat$density_est[[3]] = var(t(rbind(bs_est$density_preds[[1]],bs_est$density_preds[[2]])))/length(Y)
    var_mat[[e]][[3]] = 
      predict(main_est$fits[[e]][[1]], rbind(pred_dfs[[1]],pred_dfs[[2]]), type = "lpmatrix") %*% 
      main_est$fits[[e]][[1]]$Vp %*%
      t(predict(main_est$fits[[e]][[1]], rbind(pred_dfs[[1]],pred_dfs[[2]]), type = "lpmatrix")) + 
      predict(main_est$fits[[e]][[2]], rbind(pred_dfs[[1]],pred_dfs[[2]]), type = "lpmatrix") %*% 
      main_est$fits[[e]][[2]]$Vp %*%
      t(predict(main_est$fits[[e]][[2]], rbind(pred_dfs[[1]],pred_dfs[[2]]), type = "lpmatrix"))
    Deriv_wrt_f = as.numeric( (1/(sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]])))^2) * 
                                ( sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]])) * c(main_est$tau[[e]]$pw[[1]],main_est$tau[[e]]$pw[[2]]) - 
                                    sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]]) * c(main_est$tau[[e]]$pw[[1]],main_est$tau[[e]]$pw[[2]])) ) )
    Deriv_wrt_tau = c(main_est$density_est[[1]],main_est$density_est[[2]])
    tau_se[[e]]$agg$delta[[3]] = sqrt(t(Deriv_wrt_tau) %*% var_mat[[e]][[3]] %*% Deriv_wrt_tau + 
                                        t(Deriv_wrt_f) %*% var_mat$density_est[[3]] %*% Deriv_wrt_f)
    if (bc){
      var_mat_bc$density_est[[3]] = var_mat$density_est[[3]]
      var_mat_bc[[e]][[3]] = 
        predict(main_est$fits_bc[[e]][[1]], rbind(pred_dfs[[1]],pred_dfs[[2]]), type = "lpmatrix") %*% 
        main_est$fits_bc[[e]][[1]]$Vp %*%
        t(predict(main_est$fits_bc[[e]][[1]], rbind(pred_dfs[[1]],pred_dfs[[2]]), type = "lpmatrix")) + 
        predict(main_est$fits_bc[[e]][[2]], rbind(pred_dfs[[1]],pred_dfs[[2]]), type = "lpmatrix") %*% 
        main_est$fits_bc[[e]][[2]]$Vp %*%
        t(predict(main_est$fits_bc[[e]][[2]], rbind(pred_dfs[[1]],pred_dfs[[2]]), type = "lpmatrix"))
      Deriv_wrt_f_bc = as.numeric( (1/(sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]])))^2) * 
                                     ( sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]])) * c(main_est$tau_bc[[e]]$pw[[1]],main_est$tau_bc[[e]]$pw[[2]]) - 
                                         sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]]) * c(main_est$tau_bc[[e]]$pw[[1]],main_est$tau_bc[[e]]$pw[[2]])) ) )
      tau_se_bc[[e]]$agg$delta[[3]] = sqrt(t(Deriv_wrt_tau) %*% var_mat_bc[[e]][[3]] %*% Deriv_wrt_tau + 
                                             t(Deriv_wrt_f_bc) %*% var_mat_bc$density_est[[3]] %*% Deriv_wrt_f_bc)
    }
    if (undersmooth){
      var_mat_us$density_est[[3]] = var_mat$density_est[[3]]
      var_mat_us[[e]][[3]] = 
        predict(main_est$fits_us[[e]][[1]], rbind(pred_dfs[[1]],pred_dfs[[2]]), type = "lpmatrix") %*% 
        main_est$fits_us[[e]][[1]]$Vp %*%
        t(predict(main_est$fits_us[[e]][[1]], rbind(pred_dfs[[1]],pred_dfs[[2]]), type = "lpmatrix")) + 
        predict(main_est$fits_us[[e]][[2]], rbind(pred_dfs[[1]],pred_dfs[[2]]), type = "lpmatrix") %*% 
        main_est$fits_us[[e]][[2]]$Vp %*%
        t(predict(main_est$fits_us[[e]][[2]], rbind(pred_dfs[[1]],pred_dfs[[2]]), type = "lpmatrix"))
      Deriv_wrt_f_us = as.numeric( (1/(sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]])))^2) * 
                                     ( sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]])) * c(main_est$tau_us[[e]]$pw[[1]],main_est$tau_us[[e]]$pw[[2]]) - 
                                         sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]]) * c(main_est$tau_us[[e]]$pw[[1]],main_est$tau_us[[e]]$pw[[2]])) ) )
      tau_se_us[[e]]$agg$delta[[3]] = sqrt(t(Deriv_wrt_tau) %*% var_mat_us[[e]][[3]] %*% Deriv_wrt_tau + 
                                             t(Deriv_wrt_f_us) %*% var_mat_us$density_est[[3]] %*% Deriv_wrt_f_us)
    }
    if (mrk | kink){
      var_mat_mrk$density_est[[3]] = var_mat$density_est[[3]]
      var_mat_mrk[[e]][[3]] = 
        rbind(main_est$preds_mrk[[e]][[1]][[1]]$basis_partial_X1, main_est$preds_mrk[[e]][[2]][[1]]$basis_partial_X2) %*% 
        main_est$fits[[e]][[1]]$Vp %*%
        t(rbind(main_est$preds_mrk[[e]][[1]][[1]]$basis_partial_X1, main_est$preds_mrk[[e]][[2]][[1]]$basis_partial_X2)) + 
        rbind(main_est$preds_mrk[[e]][[1]][[2]]$basis_partial_X1, main_est$preds_mrk[[e]][[2]][[2]]$basis_partial_X2) %*% 
        main_est$fits[[e]][[2]]$Vp %*%
        t(rbind(main_est$preds_mrk[[e]][[1]][[2]]$basis_partial_X1, main_est$preds_mrk[[e]][[2]][[2]]$basis_partial_X2))
      Deriv_wrt_f_mrk = as.numeric( (1/(sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]])))^2) * 
                                      ( sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]])) * c(main_est$tau_mrk[[e]]$pw[[1]],main_est$tau_mrk[[e]]$pw[[2]]) - 
                                          sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]]) * c(main_est$tau_mrk[[e]]$pw[[1]],main_est$tau_mrk[[e]]$pw[[2]])) ) )
      tau_se_mrk[[e]]$agg$delta[[3]] = sqrt(t(Deriv_wrt_tau) %*% var_mat_mrk[[e]][[3]] %*% Deriv_wrt_tau + 
                                              t(Deriv_wrt_f_mrk) %*% var_mat_mrk$density_est[[3]] %*% Deriv_wrt_f_mrk)
      if (bc){
        var_mat_mrk_bc$density_est[[3]] = var_mat$density_est[[3]]
        var_mat_mrk_bc[[e]][[3]] = 
          rbind(main_est$preds_mrk_bc[[e]][[1]][[1]]$basis_partial_X1, main_est$preds_mrk_bc[[e]][[2]][[1]]$basis_partial_X2) %*% 
          main_est$fits_bc[[e]][[1]]$Vp %*%
          t(rbind(main_est$preds_mrk_bc[[e]][[1]][[1]]$basis_partial_X1, main_est$preds_mrk_bc[[e]][[2]][[1]]$basis_partial_X2)) + 
          rbind(main_est$preds_mrk_bc[[e]][[1]][[2]]$basis_partial_X1, main_est$preds_mrk_bc[[e]][[2]][[2]]$basis_partial_X2) %*% 
          main_est$fits_bc[[e]][[2]]$Vp %*%
          t(rbind(main_est$preds_mrk_bc[[e]][[1]][[2]]$basis_partial_X1, main_est$preds_mrk_bc[[e]][[2]][[2]]$basis_partial_X2))
        Deriv_wrt_f_mrk_bc = as.numeric( (1/(sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]])))^2) * 
                                           ( sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]])) * c(main_est$tau_mrk_bc[[e]]$pw[[1]],main_est$tau_mrk_bc[[e]]$pw[[2]]) - 
                                               sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]]) * c(main_est$tau_mrk_bc[[e]]$pw[[1]],main_est$tau_mrk_bc[[e]]$pw[[2]])) ) )
        tau_se_mrk_bc[[e]]$agg$delta[[3]] = sqrt(t(Deriv_wrt_tau) %*% var_mat_mrk_bc[[e]][[3]] %*% Deriv_wrt_tau + 
                                                   t(Deriv_wrt_f_mrk_bc) %*% var_mat_mrk_bc$density_est[[3]] %*% Deriv_wrt_f_mrk_bc)
      }
      if (undersmooth){
        var_mat_mrk_us$density_est[[3]] = var_mat$density_est[[3]]
        var_mat_mrk_us[[e]][[3]] = 
          rbind(main_est$preds_mrk_us[[e]][[1]][[1]]$basis_partial_X1, main_est$preds_mrk_us[[e]][[2]][[1]]$basis_partial_X2) %*% 
          main_est$fits_us[[e]][[1]]$Vp %*%
          t(rbind(main_est$preds_mrk_us[[e]][[1]][[1]]$basis_partial_X1, main_est$preds_mrk_us[[e]][[2]][[1]]$basis_partial_X2)) + 
          rbind(main_est$preds_mrk_us[[e]][[1]][[2]]$basis_partial_X1, main_est$preds_mrk_us[[e]][[2]][[2]]$basis_partial_X2) %*% 
          main_est$fits_us[[e]][[2]]$Vp %*%
          t(rbind(main_est$preds_mrk_us[[e]][[1]][[2]]$basis_partial_X1, main_est$preds_mrk_us[[e]][[2]][[2]]$basis_partial_X2))
        Deriv_wrt_f_mrk_us = as.numeric( (1/(sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]])))^2) * 
                                           ( sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]])) * c(main_est$tau_mrk_us[[e]]$pw[[1]],main_est$tau_mrk_us[[e]]$pw[[2]]) - 
                                               sum(c(main_est$density_preds[[1]],main_est$density_preds[[2]]) * c(main_est$tau_mrk_us[[e]]$pw[[1]],main_est$tau_mrk_us[[e]]$pw[[2]])) ) )
        tau_se_mrk_us[[e]]$agg$delta[[3]] = sqrt(t(Deriv_wrt_tau) %*% var_mat_mrk_us[[e]][[3]] %*% Deriv_wrt_tau + 
                                                   t(Deriv_wrt_f_mrk_us) %*% var_mat_mrk_us$density_est[[3]] %*% Deriv_wrt_f_mrk_us)
      }
    }
    
    ### Compute simultaneous confidence band for CATE estimates
    if (fn_ci){
      for (f in 1:2){
        max_vec[[e]][[f]] = rep(NA,reps_fn_ci)
        rho_mat = var_mat[[e]][[f]]
        for (i in 1:nrow(rho_mat)){
          rho_mat[i,] = rho_mat[i,]/sqrt(var_mat[[e]][[f]][i,i])
          rho_mat[,i] = rho_mat[,i]/sqrt(var_mat[[e]][[f]][i,i])
        }
        for (s in 1:reps_fn_ci){
          set.seed(s)
          max_vec[[e]][[f]][s] = max(abs(mvrnorm(n=nrow(pred_dfs[[f]]), mu=rep(0,nrow(pred_dfs[[f]])), Sigma = rho_mat)))
        }
        fn_ci_crit[[e]][[f]] = quantile(max_vec[[e]][[f]], probs = 1-signif_lvl)
        if (bc){
          max_vec_bc[[e]][[f]] = rep(NA,reps_fn_ci)
          rho_mat = var_mat_bc[[e]][[f]]
          for (i in 1:nrow(rho_mat)){
            rho_mat[i,] = rho_mat[i,]/sqrt(var_mat_bc[[e]][[f]][i,i])
            rho_mat[,i] = rho_mat[,i]/sqrt(var_mat_bc[[e]][[f]][i,i])
          }
          for (s in 1:reps_fn_ci){
            set.seed(s)
            max_vec_bc[[e]][[f]][s] = max(abs(mvrnorm(n=nrow(pred_dfs[[f]]), mu=rep(0,nrow(pred_dfs[[f]])), Sigma = rho_mat)))
          }
          fn_ci_crit_bc[[e]][[f]] = quantile(max_vec_bc[[e]][[f]], probs = 1-signif_lvl)
        }
        if (undersmooth){
          max_vec_us[[e]][[f]] = rep(NA,reps_fn_ci)
          rho_mat = var_mat_us[[e]][[f]]
          for (i in 1:nrow(rho_mat)){
            rho_mat[i,] = rho_mat[i,]/sqrt(var_mat_us[[e]][[f]][i,i])
            rho_mat[,i] = rho_mat[,i]/sqrt(var_mat_us[[e]][[f]][i,i])
          }
          for (s in 1:reps_fn_ci){
            set.seed(s)
            max_vec_us[[e]][[f]][s] = max(abs(mvrnorm(n=nrow(pred_dfs[[f]]), mu=rep(0,nrow(pred_dfs[[f]])), Sigma = rho_mat)))
          }
          fn_ci_crit_us[[e]][[f]] = quantile(max_vec_us[[e]][[f]], probs = 1-signif_lvl)
        }
      }
      fn_ci_crit[[e]][3] = quantile(c(max_vec[[e]][[1]],max_vec[[e]][[2]]), probs = 1-signif_lvl)
      if (bc){
        fn_ci_crit_bc[[e]][3] = quantile(c(max_vec_bc[[e]][[1]],max_vec_bc[[e]][[2]]), probs = 1-signif_lvl)
      }
      if (undersmooth){
        fn_ci_crit_us[[e]][3] = quantile(c(max_vec_us[[e]][[1]],max_vec_us[[e]][[2]]), probs = 1-signif_lvl)
      }
    }
  }
  
  for (f in 1:3){
    if (mrk & sharp_mrk){
      if (f != 3){
        tau_se_mrk[[3]]$pw[[f]] = tau_se_mrk[[1]]$pw[[f]]/abs(main_est$tau_mrk[[2]]$pw[[f]])
      }
      tau_se_mrk[[3]]$agg$delta[[f]] = tau_se_mrk[[1]]$agg$delta[[f]]/abs(main_est$tau_mrk[[2]]$agg[[f]])
      if (bc){
        if (f != 3){
          tau_se_mrk_bc[[3]]$pw[[f]] = tau_se_mrk_bc[[1]]$pw[[f]]/abs(main_est$tau_mrk_bc[[2]]$pw[[f]])
        }
        tau_se_mrk_bc[[3]]$agg$delta[[f]] = tau_se_mrk_bc[[1]]$agg$delta[[f]]/abs(main_est$tau_mrk_bc[[2]]$agg[[f]])
      }
      if (undersmooth){
        if (f != 3){
          tau_se_mrk_us[[3]]$pw[[f]] = tau_se_mrk_us[[1]]$pw[[f]]/abs(main_est$tau_mrk_us[[2]]$pw[[f]])
        }
        tau_se_mrk_us[[3]]$agg$delta[[f]] = tau_se_mrk_us[[1]]$agg$delta[[f]]/abs(main_est$tau_mrk_us[[2]]$agg[[f]])
      }
    }
    if (sharp_mrd & fuzzy){
      if (f != 3){
        tau_se[[3]]$pw[[f]] = tau_se[[1]]$pw[[f]]/abs(main_est$tau[[2]]$pw[[f]])
      }
      tau_se[[3]]$agg$delta[[f]] = tau_se[[1]]$agg$delta[[f]]/abs(main_est$tau[[2]]$agg[[f]])
      if (bc){
        if (f != 3){
          tau_se_bc[[3]]$pw[[f]] = tau_se_bc[[1]]$pw[[f]]/abs(main_est$tau_bc[[2]]$pw[[f]])
        }
        tau_se_bc[[3]]$agg$delta[[f]] = tau_se_bc[[1]]$agg$delta[[f]]/abs(main_est$tau_bc[[2]]$agg[[f]])
      }
      if (undersmooth){
        if (f != 3){
          tau_se_us[[3]]$pw[[f]] = tau_se_us[[1]]$pw[[f]]/abs(main_est$tau_us[[2]]$pw[[f]])
        }
        tau_se_us[[3]]$agg$delta[[f]] = tau_se_us[[1]]$agg$delta[[f]]/abs(main_est$tau_us[[2]]$agg[[f]])
      }
    }
  }
  
  # Compute bootstrap standard errors
  if (bs){
    for (e in 1:E){
      tau_se[[e]]$pw$bs = list(F1 = list(), F2 = list())
      tau_se[[e]]$agg$bs = list(F1 = list(), F2 = list(), all = list())
      if (bc){
        tau_se_bc[[e]]$pw$bs = tau_se[[e]]$pw$bs
        tau_se_bc[[e]]$agg$bs = tau_se[[e]]$agg$bs
      }
      if (undersmooth){
        tau_se_us[[e]]$pw$bs = tau_se[[e]]$pw$bs
        tau_se_us[[e]]$agg$bs = tau_se[[e]]$agg$bs
      }
      for (f in 1:2){
        tau_se[[e]]$pw$bs[[f]] = apply(bs_est[[e]]$pw[[f]], 1, sd)
        tau_se[[e]]$agg$bs[[f]] = sd(bs_est[[e]]$agg[[f]])
        if (bc){
          tau_se_bc[[e]]$pw$bs[[f]] = apply(bs_est_bc[[e]]$pw[[f]], 1, sd)
          tau_se_bc[[e]]$agg$bs[[f]] = sd(bs_est_bc[[e]]$agg[[f]])
        }
        if (undersmooth){
          tau_se_us[[e]]$pw$bs[[f]] = apply(bs_est_us[[e]]$pw[[f]], 1, sd)
          tau_se_us[[e]]$agg$bs[[f]] = sd(bs_est_us[[e]]$agg[[f]])
        }
      }
      tau_se[[e]]$agg$bs[[3]] = sd(bs_est[[e]]$agg[[3]])
      if (bc){
        tau_se_bc[[e]]$agg$bs[[3]] = sd(bs_est_bc[[e]]$agg[[3]])
      }
      if (undersmooth){
        tau_se_us[[e]]$agg$bs[[3]] = sd(bs_est_us[[e]]$agg[[3]])
      }
      if (mrk | kink){
        tau_se_mrk[[e]]$pw$bs = list(F1 = list(), F2 = list())
        tau_se_mrk[[e]]$agg$bs = list(F1 = list(), F2 = list(), all = list())
        if (bc){
          tau_se_mrk_bc[[e]]$pw$bs = tau_se_mrk[[e]]$pw$bs
          tau_se_mrk_bc[[e]]$agg$bs = tau_se_mrk[[e]]$agg$bs
        }
        if (undersmooth){
          tau_se_mrk_us[[e]]$pw$bs = tau_se_mrk[[e]]$pw$bs
          tau_se_mrk_us[[e]]$agg$bs = tau_se_mrk[[e]]$agg$bs
        }
        for (f in 1:2){
          tau_se_mrk[[e]]$pw$bs[[f]] = apply(bs_est_mrk[[e]]$pw[[f]], 1, sd)
          tau_se_mrk[[e]]$agg$bs[[f]] = sd(bs_est_mrk[[e]]$agg[[f]])
          if (bc){
            tau_se_mrk_bc[[e]]$pw$bs[[f]] = apply(bs_est_mrk_bc[[e]]$pw[[f]], 1, sd)
            tau_se_mrk_bc[[e]]$agg$bs[[f]] = sd(bs_est_mrk_bc[[e]]$agg[[f]])
          }
          if (undersmooth){
            tau_se_mrk_us[[e]]$pw$bs[[f]] = apply(bs_est_mrk_us[[e]]$pw[[f]], 1, sd)
            tau_se_mrk_us[[e]]$agg$bs[[f]] = sd(bs_est_mrk_us[[e]]$agg[[f]])
          }
        }
        tau_se_mrk[[e]]$agg$bs[[3]] = sd(bs_est_mrk[[e]]$agg[[3]])
        if (bc){
          tau_se_mrk_bc[[e]]$agg$bs[[3]] = sd(bs_est_mrk_bc[[e]]$agg[[3]])
        }
        if (undersmooth){
          tau_se_mrk_us[[e]]$agg$bs[[3]] = sd(bs_est_mrk_us[[e]]$agg[[3]])
        }
      }
    }
    results = list(tau = main_est$tau, tau_se = tau_se, tau_mrk = main_est$tau_mrk, bs_est = bs_est)
  }else{
    results = list(tau = main_est$tau, tau_se = tau_se, tau_mrk = main_est$tau_mrk)
  }
  if (fn_ci){
    results$fn_cov$fn_ci_crit = fn_ci_crit
  }
  if (mrk | kink){
    results$tau_se_mrk = tau_se_mrk
  }
  
  if (bc){
    results$tau_bc = main_est$tau_bc
    results$tau_se_bc = tau_se_bc
    if (bs){
      results$bs_est_bc = bs_est_bc
    }
    if (fn_ci){
      results$fn_cov$fn_ci_crit_bc = fn_ci_crit_bc
    }
    if (mrk | kink){
      results$tau_mrk_bc = main_est$tau_mrk_bc
      results$tau_se_mrk_bc = tau_se_mrk_bc
    }
  }
  
  if (undersmooth){
    results$tau_us = main_est$tau_us
    results$tau_se_us = tau_se_us
    if (bs){
      results$bs_est_us = bs_est_us
    }
    if (fn_ci){
      results$fn_cov$fn_ci_crit_us = fn_ci_crit_us
    }
    if (mrk | kink){
      results$tau_mrk_us = main_est$tau_mrk_us
      results$tau_se_mrk_us = tau_se_mrk_us
    }
  }
  
  if (save_fits){
    results$fits = main_est$fits
    if (bc){
      results$fits_bc = main_est$fits_bc
    }
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
