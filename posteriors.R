libraries <- c(
  "plyr", "ggplot2", "tidyverse", "data.table", "pracma",
  "cubature", "Rmisc", "gridExtra", "RColorBrewer", 'splines'
)
for (x in libraries) {
  library(x, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
}
#' Logit function to generate Pr_t function
#' @param t is time, in days, since infection
#' @param parameters set of model parameters as evaluated in the glm function
#' @examples
#' func_logit_cubic(100, c(1, 0.5, 0.4, .003))
#' func_logit_cubic(500, c(1, -0.5, -0.4, -.003))
func_logit_cubic <- function(t, parameters) {
  1 / (1 + exp(-(parameters[1] + parameters[2] * t + parameters[3] * t^2 + parameters[4] * t^3
  )))
}
#' cloglog function to generate Pr_t function
#' @param t is time, in days, since infection
#' @param parameters set of model parameters as evaluated in the glm function
#' @examples
#' func_cloglog_cubic(100, c(1, 0.5, 0.4, .003))
#' func_cloglog_cubic(500, c(1, -0.5, -0.4, -.003))
func_cloglog_cubic <- function(t, parameters) {
  1 - exp(-exp(parameters[1] + parameters[2] * t + parameters[3] * t^2 + parameters[4] * t^3))
}

#' Pr_t function used to generate Pr_t values using a logit link function.
#' this function receives a dataset, ODn threshold or ODn and interval between last negative to first positive.
#' @param data_set Patient dataset to be evaluated
#' @param ODn_step ODn threshold to be evaluated
#' @param t_since_inf interval size between last negative and first positive. This variable is a vector
#' @param parameters set of model parameters as evaluated in the glm function
pr_t_fun_logit_cubic <- function(data_set, ODn_step, t_since_inf) {
  pr_t <- c(0)
  dat <- data_set %>%
    mutate(eddi_1 = eddi, eddi_2 = eddi^2, eddi_3 = eddi^3)
  counter <- 0
  for (i in 1:length(ODn_step)) {
    counter <- counter + 1
    dat <- dat %>%
      mutate(recency = ifelse(ODn <= ODn_step[i] & viral_load > 1000 & eddi <= 1000, 1, 0))
    
    model <- suppressWarnings(glm2::glm2(recency ~ 1 + eddi_1 + eddi_2 + eddi_3, #
                                         family = stats::binomial(link = "logit"),
                                         data = dat
    ))
    pr_t_data <- rep(0, length(t_since_inf))
    for (j in 1:length(pr_t_data)) {
      pr_t_data[j] <- func_logit_cubic(parameters = c(
        model$coefficients[[1]], model$coefficients[[2]],
        model$coefficients[[3]], model$coefficients[[4]]
      ), t = j)
    }
    pr_t <- c(pr_t, pr_t_data)
  }
  return(pr_t[-1])
}

#' Pr_t function used to generate Pr_t values using a cloglog link function.
#' this function receives a dataset, ODn threshold or ODn and interval between last negative to first positive.
#' @param data_set Patient dataset to be evaluated
#' @param ODn_step ODn threshold to be evaluated
#' @param t_since_inf interval size between last negative and first positive. This variable is a vector
#' @param parameters set of model parameters as evaluated in the glm function
pr_t_fun_loglog_cubic <- function(data_set, ODn_step, t_since_inf) {
  pr_t <- c(0)
  dat <- data_set %>%
    mutate(eddi_1 = eddi, eddi_2 = eddi^2, eddi_3 = eddi^3)
  counter <- 0
  for (i in 1:length(ODn_step)) {
    counter <- counter + 1
    dat <- dat %>%
      mutate(recency = ifelse(ODn <= ODn_step[i] & viral_load > 1000 & eddi <= 1000, 1, 0))
    
    model <- suppressWarnings(glm2::glm2(recency ~ 1 + eddi_1 + eddi_2 + eddi_3, #
                                         family = stats::binomial(link = "cloglog"),
                                         data = dat
    ))
    
    pr_t_data <- rep(0, length(t_since_inf))
    for (j in 1:length(pr_t_data)) {
      pr_t_data[j] <- func_cloglog_cubic(parameters = c(
        model$coefficients[[1]], model$coefficients[[2]],
        model$coefficients[[3]], model$coefficients[[4]]
      ), t = j)
    }
    pr_t <- c(pr_t, pr_t_data)
  }
  return(pr_t[-1])
}

#' Generates a vector of numbers on either side of and closest to the
#' target ODn, which is used to generate Likelihood
#' @param value_x Patient ODn as measurement, after the first HIV positive test.
#' @param vec_y A vector of ODns from 0.01 to 4 step size 0.01
#' @param to_select number to be selected from of ODns closest of the target ODn, to be selected.
nearest_numbers <- function(value_x, vec_y, to_select) {
  # browser()
  y <- sort(vec_y)
  value_diff <- abs(value_x - y)
  dat <- (data.frame(y, value_diff) %>%
            arrange(value_diff) # %>%
          # filter(y!= value_x)
  )[1:(to_select + 1), 1]
  dat_1 <- ifelse(dat %in% value_x, dat, c(dat[1:to_select], value_x))
  return(dat_1)
}

#' Generates the glm model parameters used in the computation of likelihood of the target ODn
#' given an inter-test interval target ODn
#' @param dat a dataset of Pr_t values by ODn threshold and per day over 1000 daystime that's brought in as a matrix.
#' @param value_x Patient ODn as measurement, after the first HIV positive test.
#' @param to_select number to be selected from of ODns closest of the target ODn, to be selected.
#' @param t_since_ln duration of inter-test interval.
#'
likelihood_param_quad_function <- function(dat, target_ODn, target_ODnstepsize, t_since_ln) {
  vec_x <- c(nearest_numbers(value_x = target_ODn, vec_y = seq(0.01, 5, target_ODnstepsize), to_select = 6))
  dat <- subset(as.data.frame(dat), threshold %in% vec_x)
  pr_t_slope_data <- data.frame(
    t_var = NA,
    intercept = NA,
    linear_term = NA,
    quad_term = NA
  )
  counter <- 0
  for (i in 1:length(t_since_ln)) {
    counter <- counter + 1
    m1 <- suppressWarnings(glm2::glm2(pr_t ~ 1 + threshold + I(threshold^2), #
                                      family = stats::gaussian(link = "identity"),
                                      data = subset(dat, vec_time == t_since_ln[i])
    ))
    
    pr_t_slope_data[counter, ] <- c(
      t_var = vec_time[i],
      intercept = m1$coefficients[[1]],
      linear_term = m1$coefficients[[2]],
      quad_term = m1$coefficients[[3]]
    )
  }
  return(pr_t_slope_data)
}

#' Generates the likelihood of the target ODn given an inter-test interval.
#' We differentiate a quadratic polynomial and evaluate it, per time point in the inter-test interval, using
#' parameters generated in likelihood_param_quad_function.
#' @param param_datset dataset of intercept (not really important as it doesn't get used),
#' linear, and quadratic terms that are used in evaluating the quad function.
#' @param ODn is the target ODn
#' @param t_since_ln vector of time points in the inter-test interval
#'
likelihood_fun <- function(param_datset, ODn, t_since_ln) {
  l <- rep(0, length(t_since_ln))
  ff <- expression(1 + param1 * t + param2 * t^2)
  f <- D(ff, "t")
  counter <- 0
  for (j in 1:length(t_since_ln)) {
    counter <- counter + 1
    param1 <- as.numeric(param_datset$linear_term[[j]])
    param2 <- as.numeric(param_datset$quad_term[[j]])
    t <- ODn
    l[counter] <- eval(f)
  }
  likelihood <- data.frame(l) %>%
    mutate(
      time_t = t_since_ln,
      bigL = l / (trapz(time_t, l))
    )
  return(likelihood)
}

# tic()
# pr_t_logit_cubic <- data.frame(pr_t = pr_t_fun_logit_cubic(data_set = data_generate_pr_t, ODn_step = ODn_step, t_since_inf = vec_time)) %>%
#   mutate(
#     threshold = rep(ODn_step, each = length(vec_time)),
#     vec_time = rep(vec_time, times = length(ODn_step))
#   )
# toc()
#
# tic()
# pr_t_loglog_cubic <- data.frame(pr_t = pr_t_fun_loglog_cubic(data_set = data_generate_pr_t, ODn_step = ODn_step, t_since_inf = vec_time)) %>%
#   mutate(
#     threshold = rep(ODn_step, each = length(vec_time)),
#     vec_time = rep(vec_time, times = length(ODn_step))
#   )
# toc()

# saveRDS(pr_t_logit_cubic, "pr_t_logit_evaluations.rds")
# saveRDS(pr_t_loglog_cubic, "pr_t_loglog_evaluations.rds")

pr_t_logit_cubic <- readRDS("pr_t_logit_evaluations.rds")
# pr_t_loglog_cubic <- readRDS("pr_t_loglog_evaluations.rds")