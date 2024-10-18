# start by simulating trial data with a time varying foi

library(tidyverse)

# weibull fn
weibull <- function(time, alpha, beta) {
  pow <- -(time/beta)^alpha
  y <- exp(pow)
  return(y)
}

# generate a random walk
# mu = drift
RW <- function(len, x0, mu = 0, variance) {
  z<-cumsum(rnorm(n=len, mean=0, 
                  sd=sqrt(variance)))
  t<-1:len
  x<-x0+t*mu+z
  return(x)
}

quantile_95 <- function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
}

# function to simualate a trial arm
simulate_trial_arm <- function(N_patients, max_time, 
                               drug_arm, 
                               alpha, beta, lambda_vec) {
  times <- seq(1, max_time, by = 1)
  if(drug_arm == 0) {
    prob_susceptible <- rep(1, length(times))
  } else {
    prob_susceptible <- 1 - weibull(time = times,
                                    alpha = alpha,
                                    beta = beta)
  }
  
  lambda <- rep(lambda_vec, each = 7)[1:max_time]
  
  N_remaining <- N_patients
  df <- data.frame(days_since_dose = times,
                   n_inf =  NA,
                   cum_inf =  NA,
                   n_sus =  NA,
                   arm = drug_arm)
  
  for(i in 1:nrow(df)) {
    df$n_sus[i] <- N_remaining
    df$n_inf[i] <- rbinom(1, N_remaining, 1 - exp(-prob_susceptible[i] * lambda[i]))
    df$cum_inf[i] <- sum(df$n_inf[1:i])
    N_remaining <- N_patients - df$cum_inf[i]
  }
  
  return(df)
}


# set up the MCMC
# do we need lambdas in here too?
df_params <- define_params(name = "mu", min = 0, max = 1, # RW mean
                           name = "sigma", min = 0, max = 1, # RW var
                           name = "lambda_1", min = 0, max = 1,
                           name = "lambda_2", min = 0, max = 1,
                           name = "lambda_3", min = 0, max = 1,
                           name = "lambda_4", min = 0, max = 1,
                           name = "lambda_5", min = 0, max = 1,
                           name = "alpha", min = 1, max = 10,
                           name = "beta", min = 0, max = 100)
print(df_params)

r_loglike <- function(params, data, misc) {
  
  # extract parameter values
  lambda_1 <- params["lambda_1"]
  lambda_2 <- params["lambda_2"]
  lambda_3 <- params["lambda_3"]
  lambda_4 <- params["lambda_4"]
  lambda_5 <- params["lambda_5"]
  alpha_spaq <- params["alpha_spaq"]
  beta_spaq <- params["beta_spaq"]
  alpha_dhapq <- params["alpha_dhapq"]
  beta_dhapq <- params["beta_dhapq"]
  min <- 1e-10 # to prevent zeroing out the likelihood when infections with 100% weibull PE
  
  lambda_spaq <- lambda * ((1-min)*(1 - weibull(time = data$t_spaq, 
                                                alpha = alpha_spaq,
                                                beta = beta_spaq)) + min)
  lambda_dhapq <- lambda * ((1-min)*(1 - weibull(time = data$t_dhapq, 
                                                 alpha = alpha_dhapq,
                                                 beta = beta_dhapq)) + min)
  
  ret <- 
    sum(dbinom(x = data$control_n_inf,
               size = data$control_n_sus,
               prob = 1 - exp(-lambda),
               log = TRUE)) +
    sum(dbinom(x = data$spaq_n_inf,
               size = data$spaq_n_sus,
               prob = 1 - exp(-lambda_spaq),
               log = TRUE)) +
    sum(dbinom(x = data$dhapq_n_inf,
               size = data$dhapq_n_sus,
               prob = 1 - exp(-lambda_dhapq),
               log = TRUE))
  
  return(ret)
}

r_logprior <- function(params, misc) {
  
  # extract parameter values
  lambda_1 <- params["lambda_1"]
  lambda_2 <- params["lambda_2"]
  lambda_3 <- params["lambda_3"]
  lambda_4 <- params["lambda_4"]
  lambda_5 <- params["lambda_5"]
  mu <- params["mu"]
  sigma2 <- params["sigma2"]
  alpha_spaq <- params["alpha_spaq"]
  beta_spaq <- params["beta_spaq"]
  alpha_dhapq <- params["alpha_dhapq"]
  beta_dhapq <- params["beta_dhapq"]
  
  ## gamma based upon prior sampling
  ret <- dnorm(lambda_1, mu, sqrt(sigma2)) +
    dnorm(lambda_2, lambda_1, sqrt(sigma2)) +
    dnorm(lambda_3, lambda_2, sqrt(sigma2)) +
    dnorm(lambda_4, lambda_3, sqrt(sigma2)) +
    dnorm(lambda_5, lambda_4, sqrt(sigma2)) +
    dgamma(alpha_spaq, 4, 1, log = TRUE) +
    dgamma(alpha_dhapq, 4, 1, log = TRUE) +
    dgamma(beta_spaq, 10, 0.3, log = TRUE) +
    dgamma(beta_dhapq, 10, 0.3, log = TRUE)
  
  
  return(ret)
}

# simulate trials to test
N <- 3e3
time <- 30
t_vec <- 1:time
mu <- 1e-2
sigma2 <- 1e-5
n_trials <- 100

# generate a simulated set of parameters
set.seed(123)
real_params <- data.frame(alpha_spaq = rnorm(n_trials, 4, 1),
                          beta_spaq = rnorm(n_trials, 40, 0.2),
                          alpha_dhapq = rnorm(n_trials, 4, 1),
                          beta_dhapq = rnorm(n_trials, 40, 0.2))
for(i in 1:nrow(real_params)) {
  # need to draw each simulation foi independently
  lambda <- RW(len = 5, x0 = mu, variance = sigma2)
  real_params$lambda_1[i] <- lambda[1]
  real_params$lambda_2[i] <- lambda[2]
  real_params$lambda_3[i] <- lambda[3]
  real_params$lambda_4[i] <- lambda[4]
  real_params$lambda_5[i] <- lambda[5]
}

patients <- list(control = 500,
                 spaq = 2000,
                 dhapq = 2000)

control <- simulate_trial_arm(N_patients = patients$control,
                              max_time = time,
                              drug_arm = 0,
                              lambda_vec = c(real_params$lambda_1[i],
                                             real_params$lambda_2[i],
                                             real_params$lambda_3[i],
                                             real_params$lambda_4[i],
                                             real_params$lambda_5[i]))
spaq <- simulate_trial_arm(N_patients = patients$spaq,
                           max_time = time,
                           drug_arm = 1,
                           alpha = real_params$alpha_spaq[i],
                           beta = real_params$beta_spaq[i],
                           lambda_vec = c(real_params$lambda_1[i],
                                          real_params$lambda_2[i],
                                          real_params$lambda_3[i],
                                          real_params$lambda_4[i],
                                          real_params$lambda_5[i]))
dhapq <- simulate_trial_arm(N_patients = patients$dhapq,
                            max_time = time,
                            drug_arm = 2,
                            alpha = real_params$alpha_dhapq[i],
                            beta = real_params$beta_dhapq[i],
                            lambda_vec = c(real_params$lambda_1[i],
                                           real_params$lambda_2[i],
                                           real_params$lambda_3[i],
                                           real_params$lambda_4[i],
                                           real_params$lambda_5[i]))
out <- rbind(control, spaq, dhapq)
out$trial <- i
assign(paste0("sim_",i), out)
assign(paste0("sim_list_", i), list(control_n_sus = out$n_sus[out$arm == 0],
                                    control_n_inf = out$n_inf[out$arm == 0],
                                    t_control = out$days_since_dose[out$arm == 0],
                                    spaq_n_sus = out$n_sus[out$arm == 1],
                                    spaq_n_inf = out$n_inf[out$arm == 1],
                                    t_spaq = out$days_since_dose[out$arm == 1],
                                    dhapq_n_sus = out$n_sus[out$arm == 2],
                                    dhapq_n_inf = out$n_inf[out$arm == 2],
                                    t_dhapq = out$days_since_dose[out$arm == 2]))
# print(paste0("data created for iteration ", i, "-- start mcmc run"))
out_mcmc <- run_mcmc(data = list(control_n_sus = out$n_sus[out$arm == 0],
                                 control_n_inf = out$n_inf[out$arm == 0],
                                 t_control = out$days_since_dose[out$arm == 0],
                                 spaq_n_sus = out$n_sus[out$arm == 1],
                                 spaq_n_inf = out$n_inf[out$arm == 1],
                                 t_spaq = out$days_since_dose[out$arm == 1],
                                 dhapq_n_sus = out$n_sus[out$arm == 2],
                                 dhapq_n_inf = out$n_inf[out$arm == 2],
                                 t_dhapq = out$days_since_dose[out$arm == 2]),
                     df_params = df_params, 
                     loglike = r_loglike,
                     logprior = r_logprior,
                     burnin = 1e3,
                     samples = 10e3,
                     chains = 5,
                     pb_markdown = TRUE)
