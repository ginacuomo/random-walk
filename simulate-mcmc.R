## try and simulate data for a control arm from a random walk and then fit an MCMC returning the mean and sd

simulate_trial_arm <- function(N_patients, 
                               max_time, 
                               # drug_arm, 
                               lambda_vec) {
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
    df$n_inf[i] <- rbinom(1, N_remaining, 1 - exp(-prob_susceptible[i] * lambda))
    df$cum_inf[i] <- sum(df$n_inf[1:i])
    N_remaining <- N_patients - df$cum_inf[i]
  }
  
  return(df)
}

# loglikelihood
r_loglike <- function(params, data, misc) {
  
  # extract parameter values
  lambda_1 <- params["lambda_1"]
  lambda_2 <- params["lambda_2"]
  lambda_3 <- params["lambda_3"]
  lambda_4 <- params["lambda_4"]
  lambda_5 <- params["lambda_5"]
  lambda_6 <- params["lambda_6"]
  lambda_7 <- params["lambda_7"]
  lambda_8 <- params["lambda_8"]
  lambda <- c(rep(c(lambda_1, lambda_2, lambda_3, lambda_4, 
                    lambda_5, lambda_6, lambda_7, lambda_8), 
                  each = 7))[1:max_time]
  
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

r_logprior <- function(params, misc) {
  
  # extract parameter values
  lambda_1 <- params["lambda_1"]
  lambda_2 <- params["lambda_2"]
  lambda_3 <- params["lambda_3"]
  lambda_4 <- params["lambda_4"]
  lambda_5 <- params["lambda_5"]
  lambda_6 <- params["lambda_6"]
  lambda_7 <- params["lambda_7"]
  lambda_8 <- params["lambda_8"]
  mu <- params["mu"]
  sigma2 <- params["sigma2"]
  
  # normal distribution based on a 
  ret <- 
    dnorm(lambda_1, mu, sqrt(sigma2)) + 
    dnorm(lambda_2, lambda_1, sqrt(sigma2)) + 
    dnorm(lambda_3, lambda_2, sqrt(sigma2)) + 
    dnorm(lambda_4, lambda_3, sqrt(sigma2)) + 
    dnorm(lambda_5, lambda_4, sqrt(sigma2)) + 
    dnorm(lambda_6, lambda_5, sqrt(sigma2)) + 
    dnorm(lambda_7, lambda_6, sqrt(sigma2)) + 
    dnorm(lambda_8, lambda_7, sqrt(sigma2))
  return(ret)
}

df_params <- define_params(name = "lambda_1", min = 1e-3, max = 1,
                           name = "lambda_2", min = 1e-3, max = 1,
                           name = "lambda_3", min = 1e-3, max = 1,
                           name = "lambda_4", min = 1e-3, max = 1,
                           name = "lambda_5", min = 1e-3, max = 1,
                           name = "lambda_6", min = 1e-3, max = 1,
                           name = "lambda_7", min = 1e-3, max = 1,
                           name = "lambda_8", min = 1e-3, max = 1,
                           name = "mu", min = 1e-3, max = 1,
                           name = "sigma2", min = 1e-10, max = 3)

# simulate a trial
RW <- function(N, x0, mu, variance) {
  z <- cumsum(rnorm(n = N, mean = 0, 
                    sd = sqrt(variance)))
  t <- 1:N
  x <- x0 + t*mu + z
  return(x)
}

set.seed(123)
foi <- RW(8, x0 = 0.01, mu = 0, variance = 1e-5)
control <- simulate_trial_arm(1000, 56, foi)

start <- Sys.time()
out_mcmc <- run_mcmc(data = list(control_n_sus = control$n_sus,
                                 control_n_inf = control$n_inf),
                     df_params = df_params, 
                     loglike = r_loglike,
                     logprior = r_logprior,
                     burnin = 1e3,
                     samples = 1e3,
                     chains = 5,
                     silent = TRUE)
Sys.time() - start
q95 <- out_mcmc$output %>%
  filter(phase == "sampling") %>%
  select(-c(chain, phase, iteration, logprior, loglikelihood)) %>%
  apply(2, quantile_95) %>%
  as.data.frame()

## now simulate 100 control datasets and estimate the power with which it estimates mu and sigma2
n_trials <- 1e2
real_params <- data.frame(lambda_1 = rep(0, n_trials))
set.seed(1234)
for(i in 1:nrow(real_params)) {
  mu <- rnorm(n = 1, mean = 0.03, sd = 0.008)
  sigma2 <- runif(n = 1, min = 5e-6, max = 2e-5)
  foi <- RW(8, x0 = 0.03, mu = 0, variance = 1e-5)
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
  q95 <- out_mcmc$output %>%
    filter(phase == "sampling") %>%
    select(-c(chain, phase, iteration, logprior, loglikelihood)) %>%
    apply(2, quantile_95) %>%
    as.data.frame()
  
  # wrangle and store results
  res_list[[i]] <- t(q95) %>%
    as.data.frame() %>%
    mutate(true = c(unlist(real_params[i,])),
           sim = i,
           param = names(q95))
}
Sys.time() - t0







