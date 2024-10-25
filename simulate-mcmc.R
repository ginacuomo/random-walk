# simulate-mcmc.R
#
# Author: Gina Cuomo-Dannenburg
# Date: 2024-10-23
#
# Inputs: (none)
#
# Outputs: (none)
#
# Purpose:
# Try and simulate data for a control arm from a random walk and then fit an
# MCMC returning the mean and sd.
#
# ------------------------------------------------------------------
# devtools::install_github("mrc-ide/drjacoby@v1.5.4")
library(drjacoby)
library(dplyr)
library(ggplot2)

# draw from a random walk with drift
RW <- function(N, x0, mu, variance) {
  z <- rep(0, N)
  z[1] <- abs(rnorm(1, mean = x0, sd = sqrt(variance)))
  for (i in 2:N) {
    z[i] <- abs(rnorm(1, mean = z[i-1] + mu, sd = sqrt(variance)))
  }
  return(z)
}

# simulate one arm of a trial
simulate_trial_arm <- function(N_patients, 
                               max_time, 
                               lambda_vec,
                               drug_arm = NA) {
  
  times <- seq(1, max_time, by = 1)
  # if(drug_arm == 0) {
  prob_susceptible <- rep(1, length(times))
  # } else {
  #   prob_susceptible <- 1 - weibull(time = times,
  #                                   alpha = alpha,
  #                                   beta = beta)
  # }
  
  N_remaining <- N_patients
  df <- data.frame(days_since_dose = times,
                   n_inf =  NA,
                   cum_inf =  NA,
                   n_sus =  NA,
                   arm = drug_arm)
  
  lambda <- rep(lambda_vec, each = 7)[1:max_time]
  
  for(i in 1:nrow(df)) {
    df$n_sus[i] <- N_remaining
    df$n_inf[i] <- rbinom(1, N_remaining, 1 - exp(-prob_susceptible[i] * lambda[i]))
    df$cum_inf[i] <- sum(df$n_inf[1:i])
    N_remaining <- N_patients - df$cum_inf[i]
  }
  
  return(df)
}

# loglikelihood
r_loglike <- function(params, data, misc) {
  
  # extract parameter values
  params_lambda <- params[sprintf("lambda_%s", 1:8)]
  lambda <- rep(params_lambda, each = 7)[1:max_time]
  
  ret <- 
    sum(dbinom(x = data$control_n_inf,
               size = data$control_n_sus,
               prob = 1 - exp(-lambda),
               log = TRUE)) 
  
  return(ret)
}

# logprior
r_logprior <- function(params, misc) {
  
  # extract parameter values
  params_lambda <- params[sprintf("lambda_%s", 1:8)]
  x0 <- params["x0"]
  mu <- params["mu"]
  sigma2 <- params["sigma2"]
  
  # normal distribution based on a random walk
  ret <- dnorm(params_lambda[1], mean = x0, sd = sqrt(sigma2), log = TRUE)
  for (i in 2:length(params_lambda)) {
    ret <- ret + dnorm(params_lambda[i], 
                       mean = params_lambda[i-1] + mu, 
                       sd = sqrt(sigma2), log = TRUE)
  }
  
  return(ret)
}

# setup parameters data.frame
df_params <- data.frame(name = sprintf("lambda_%s", 1:8), min = 0, max = 1) |>
  bind_rows(define_params(name = "x0", min = 0, max = 1,
                          name = "mu", min = 0, max = 1,
                          name = "sigma2", min = 0, max = 1))

# simulate a trial
#set.seed(123)
foi <- RW(8, x0 = 0.05, mu = 0.01, variance = 1e-4)
control <- simulate_trial_arm(N_patients = 1000,
                              max_time = 56,
                              lambda_vec = foi)
max_time <- 56

# quick KM plot of trial
plot(control$days_since_dose, control$n_sus, type = "s", ylim = c(0, max(control$n_sus)))

# run MCMC
start <- Sys.time()
out_mcmc <- run_mcmc(data = list(control_n_sus = control$n_sus,
                                 control_n_inf = control$n_inf),
                     df_params = df_params, 
                     loglike = r_loglike,
                     logprior = r_logprior,
                     burnin = 1e3,
                     samples = 1e3,
                     chains = 5)
Sys.time() - start

# plot CrIs on lambda
lambda_names <- sprintf("lambda_%s", 1:8)
drjacoby::plot_credible(out_mcmc, show = lambda_names) +
  geom_point(aes(x = param, y = value), data.frame(param = lambda_names, value = foi), col = "red")
# get CrIs
q95 <- out_mcmc$output |>
  filter(phase == "sampling") |>
  select(-c(chain, phase, iteration, logprior, loglikelihood)) |>
  apply(2, drjacoby:::quantile_95) |>
  as.data.frame()

# POWER ANALYSIS
# ------------------------------------------------------------------------------------------------
# TODO: check the syntax switch from mu --> x0

## now simulate 100 control datasets and estimate the power with which it estimates mu and sigma2
n_trials <- 1e2
real_params <- data.frame(lambda_1 = rep(0, n_trials))
for(i in 1:nrow(real_params)) {
  x0 <- rnorm(n = 1, mean = 0.03, sd = 0.008)
  mu <- rnorm(n = 1, mean = 0.01, sd = 1e-4)
  sigma2 <- runif(n = 1, min = 5e-6, max = 2e-5)
  foi <- RW(8, x0 = 0.03, mu = mu, variance = sigma2)
  real_params$lambda_1[i] <- foi[1]
  real_params$lambda_2[i] <- foi[2]
  real_params$lambda_3[i] <- foi[3]
  real_params$lambda_4[i] <- foi[4]
  real_params$lambda_5[i] <- foi[5]
  real_params$lambda_6[i] <- foi[6]
  real_params$lambda_7[i] <- foi[7]
  real_params$lambda_8[i] <- foi[8]
  real_params$x0[i] <- x0
  real_params$mu[i] <- mu
  real_params$sigma2[i] <- sigma2
}
sum(real_params<0)
res_list <- list()
t0 <- Sys.time()
for (i in 1:n_trials) {
  message(sprintf("i = %s", i))
  
  # simulate data
  control <- simulate_trial_arm(1000, 56, lambda_vec = c(real_params$lambda_1[i],
                                                         real_params$lambda_2[i],
                                                         real_params$lambda_3[i],
                                                         real_params$lambda_4[i],
                                                         real_params$lambda_5[i],
                                                         real_params$lambda_6[i],
                                                         real_params$lambda_7[i],
                                                         real_params$lambda_8[i]))
  
  # run MCMC
  out_mcmc <- run_mcmc(data = list(control_n_sus = control$n_sus,
                                   control_n_inf = control$n_inf),
                       df_params = df_params, 
                       loglike = r_loglike,
                       logprior = r_logprior,
                       burnin = 1e3,
                       samples = 1e3,
                       chains = 5,
                       silent = TRUE)
  
  # plot_par(out_mcmc, phase = "both")
  
  # get CrIs
  q95 <- out_mcmc$output |>
    filter(phase == "sampling") |>
    select(-c(chain, phase, iteration, logprior, loglikelihood)) |>
    apply(2, drjacoby:::quantile_95) |>
    as.data.frame()
  
  # wrangle and store results
  res_list[[i]] <- t(q95) |>
    as.data.frame() |>
    mutate(true = c(unlist(real_params[i,])),
           sim = i,
           param = names(q95))
}
Sys.time() - t0

ret <- res_list |>
  bind_rows()
row.names(ret) <- NULL
ret$inside <- (ret$true > ret$Q2.5) & (ret$true < ret$Q97.5)

# count times correct - this is the power to accurately estimate parameters 
ret_count <- ret |>
  group_by(param) |>
  summarise(inside = sum(inside)) |>
  dplyr::rowwise() |>
  dplyr::filter(grepl("lambda", param))

ggplot(data = ret_count, (aes(x = param, y = inside))) + geom_point() + 
  geom_hline(yintercept = 95, lty = 2) + theme_bw() + ylim(c(90, 100))

