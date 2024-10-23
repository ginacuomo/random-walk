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
# BOB CHANGES TO CODE
# - in simulate_trial_arm, given drug_arm default value NA rather than commenting out, as this parameter is used in the function
# - in simulate_trial_arm, the expression rbinom(1, N_remaining, 1 - exp(-prob_susceptible[i] * lambda)) was missing the lambda[i] index, meaning it won't pull the correct values.
# - changed RW function so values are reflected around zero, i.e. cannot generate negative values
# - swapped pipe symbol from %>% to |> (just because I'm pedantic!)
# - simplified definition and extraction of values lambda_1, lambda_2 etc. using sprintf
# - random walk simulator (RW function) used x0 as intercept and mu as drift term. But prior used mu as intercept with no drift term. Modified params data.frame and prior function to have x0 intercept and mu drift term
# - random walk simulator used t*mu drift term, where t=1:N. But this applies drift even to the first observation, i.e. intercept becomes x0 + mu. Changed so t effectively starts at 0 meaning intercept is x0 (easier to interpret)
# - changed range of some parameters in df_params to start at zero, rather than starting at arbitrarily low values. It's fine for the range to go all the way down to zero, we don't have to worry about the MCMC reaching values of exactly zero (this will never happen)
# - THE MOST IMPORTANT CHANGE! logprior forgot to log the normal density (we've all been here)!
#
# ------------------------------------------------------------------
#devtools::install_github("mrc-ide/drjacoby@v1.5.4")
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
               log = TRUE)) #+
  # sum(dbinom(x = data$spaq_n_inf,
  #            size = data$spaq_n_sus,
  #            prob = 1 - exp(-lambda_spaq),
  #            log = TRUE)) +
  # sum(dbinom(x = data$dhapq_n_inf,
  #            size = data$dhapq_n_sus,
  #            prob = 1 - exp(-lambda_dhapq),
  #            log = TRUE))
  
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
    ret <- ret + dnorm(params_lambda[i], mean = params_lambda[i-1] + mu, sd = sqrt(sigma2), log = TRUE)
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
                     chains = 1)
Sys.time() - start

# plot CrIs on lambda
lambda_names <- sprintf("lambda_%s", 1:8)
drjacoby::plot_credible(out_mcmc, show = lambda_names) +
  geom_point(aes(x = param, y = value), data.frame(param = lambda_names, value = foi), col = "red")

# ------------------------------------------------------------------------------------------------
# BOB STOPPED EDITING

# get CrIs
q95 <- out_mcmc$output |>
  filter(phase == "sampling") |>
  select(-c(chain, phase, iteration, logprior, loglikelihood)) |>
  apply(2, drjacoby:::quantile_95) |>
  as.data.frame()

## now simulate 100 control datasets and estimate the power with which it estimates mu and sigma2
n_trials <- 1e2
real_params <- data.frame(lambda_1 = rep(0, n_trials))
set.seed(1234)
for(i in 1:nrow(real_params)) {
  mu <- rnorm(n = 1, mean = 0.03, sd = 0.008)
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
    apply(2, quantile_95) |>
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
ret$inside <- (ret$true > ret$`2.5%`) & (ret$true < ret$`97.5%`)

# count times correct
ret_count <- ret |>
  group_by(param) |>
  summarise(inside = sum(inside))

exp <- ret |> dplyr::filter(param == "sigma2") |>
  dplyr::mutate(estimate_div_true = `50%`/true) 
mean(exp$estimate_div_true)


