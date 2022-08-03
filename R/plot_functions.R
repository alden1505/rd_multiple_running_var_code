#' Subfunction for plotting
#' 
#' default_ggplot is a subfunction that takes in a ggplot dataframe and
#' a number of arguments and creates a ggplot object.
#'
#' @import locfit
#' @import mgcv
#' @import Hmisc
#' @import ggplot2
#' @param df ggplot dataframe
#' @param my.colors specify desired colors for plot
#' @param my_ylab title for y-axis
#' @param my_xlab title for x-axis
#' @param plot_ci specify whether we want to plot confidence intervals, defaults to TRUE
#' @keywords plot_subfn
#' @export

default_ggplot = function(df, my.colors, my_ylab, my_xlab, plot_ci = TRUE){
  my_plot = ggplot(data = df, aes(x = xaxis, y = treatment,
                                  group = type, shape = type,
                                  colour = type, fill = type)) +
    geom_line(aes(linetype = type), size = 1) +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    scale_color_manual(values = my.colors) +
    scale_fill_manual(values = my.colors) +
    scale_linetype_manual(values = c(1,1,0)) +
    xlim(c( min(df$x), max(df$x) )) +
    ylim( c( min(c(0, df$lower, df$treatment), na.rm = TRUE),
             max(c(0, df$upper, df$treatment), na.rm = TRUE) ) ) +
    theme_classic() +
    theme(axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black")) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=16)) + 
    labs(y = my_ylab, x = my_xlab)
  if (plot_ci){
    my_plot + geom_ribbon(aes(x = xaxis, ymin = lower, ymax = upper,
                              fill = type), linetype = 0,  alpha = 0.1)
  }
}

#' Function to create basic MRD plot
#'
#' create_plot is a function that takes in the range for a running variable, 
#' estimates of the treatment effect and standard errors along that range, and other arguments, 
#' and returns a ggplot object, using the default_ggplot subfunction.
#' This function is mainly used to generate basic plot for the MRD estimates.
#' 
#' @import locfit
#' @import mgcv
#' @import Hmisc
#' @import ggplot2
#' @param x range for running variable
#' @param est estimates of the treatment effects
#' @param se standard error of estimated treatment effects
#' @param legend title for legend
#' @param reverse reverse the running variable, defaults to FALSE
#' @param my.colors specify desired colors for plot, defaults to black
#' @param my_xlab title for x-axis
#' @param my_ylab title for y-axis
#' @keywords basic_plot
#' @export

create_plot = function(x, est, se,
                       legend = "MRD", reverse = FALSE,
                       my.colors = "black",
                       my_xlab, my_ylab){
  if (reverse){
    df = data.frame(xaxis = x,
                    treatment = -est,
                    lower = -est - qnorm(0.975)*se,
                    upper = -est + qnorm(0.975)*se,
                    type = factor (rep( legend, each = length(x) ),
                                   levels = legend ) )
  }else{
    df = data.frame(xaxis = x,
                    treatment = est,
                    lower = est - qnorm(0.975)*se,
                    upper = est + qnorm(0.975)*se,
                    type = factor (rep( legend, each = length(x) ),
                                   levels = legend ) )
  }
  my_plot = default_ggplot(df, my.colors, my_ylab, my_xlab)
  return(my_plot)
}

# This function takes in a ggplot object created by create_plot,
# and adds point estimates and standard errors from a different estimator.
add_plot = function(ggplot_obj, new_est = NA, new_se = NA,
                    new_ci_lower = NA, new_ci_upper = NA,
                    plot_ci = FALSE,
                    reverse = FALSE, constant = TRUE,
                    x = NA, est_obj = NA, rep_times = NA,
                    new_label, my.colors){
  if (constant){
    if (reverse){
      new_est = -new_est
      new_ci_lower_reverse = -new_ci_upper
      new_ci_upper_reverse = -new_ci_lower
      new_ci_lower = new_ci_lower_reverse
      new_ci_upper = new_ci_upper_reverse
    }
    if (!plot_ci){
      ggplot_df = rbind(ggplot_obj$data,
                        data.frame(xaxis = ggplot_obj$data$xaxis[ggplot_obj$data$type == "MRD"],
                                   treatment = new_est,
                                   lower = new_est - qnorm(0.975)*new_se,
                                   upper = new_est + qnorm(0.975)*new_se,
                                   type = new_label))
    }else{
      ggplot_df = rbind(ggplot_obj$data,
                        data.frame(xaxis = ggplot_obj$data$xaxis[ggplot_obj$data$type == "MRD"],
                                   treatment = NA,
                                   lower = new_ci_lower,
                                   upper = new_ci_upper,
                                   type = new_label))
    }
  }else{
    x = seq(from = min(x), to = max(x), length = rep_times*(length(x)-1))
    new_est = numeric(length(est_obj))
    new_se = numeric(length(est_obj))
    new_ci_lower = numeric(length(est_obj))
    new_ci_upper = numeric(length(est_obj))
    for (i in 1:length(est_obj)){
      new_est[i] = est_obj[[i]]$point_est
      new_se[i] = est_obj[[i]]$se
      new_ci_lower[i] = est_obj[[i]]$ci_lower
      new_ci_upper[i] = est_obj[[i]]$ci_upper
    }
    if (reverse){
      new_est = -new_est
      new_ci_lower_reverse = -new_ci_upper
      new_ci_upper_reverse = -new_ci_lower
      new_ci_lower = new_ci_lower_reverse
      new_ci_upper = new_ci_upper_reverse
    }
    if (!plot_ci){
      ggplot_df = rbind(ggplot_obj$data,
                        data.frame(xaxis = x,
                                   treatment = rep(new_est, each=rep_times),
                                   lower = rep(new_est, each=rep_times) - qnorm(0.975)*rep(new_se, each=rep_times),
                                   upper = rep(new_est, each=rep_times) + qnorm(0.975)*rep(new_se, each=rep_times),
                                   type = new_label))
    }else{
      ggplot_df = rbind(ggplot_obj$data,
                        data.frame(xaxis = x,
                                   treatment = NA,
                                   lower = rep(new_ci_lower, each=rep_times),
                                   upper = rep(new_ci_upper, each=rep_times),
                                   type = new_label))
    }
  }
  default_ggplot(df = ggplot_df,
                 my.colors = my.colors,
                 my_ylab = ggplot_obj$labels$y,
                 my_xlab = ggplot_obj$labels$x)
}
