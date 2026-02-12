# stan_code <- "
# data {
#   int<lower=1> T;
#   int<lower=1> P;
#   int<lower=1> N_obs;
# 
#   int<lower=1,upper=T> obs_t[N_obs];
#   vector[N_obs] y_obs;
# 
#   int<lower=1,upper=P> period[T];
# 
#   vector[P] beta_mu;
#   vector[P] beta_se;
# }
# parameters {
#   vector[T] y;
#   vector[P] beta;
# 
#   real<lower=0> sigma_proc;
#   real<lower=0> sigma_obs;
# }
# model {
#   beta ~ normal(beta_mu, beta_se);
# 
#   sigma_proc ~ exponential(1);
#   sigma_obs  ~ exponential(1);
# //  sigma_proc ~ student_t(5, 0, 2);
# //  sigma_obs  ~ student_t(10, 0, 2);
# //    sigma_proc ~ normal(0, 2);
# //    sigma_obs  ~ normal(0, 2);
# 
#   y[1] ~ normal(355, 8.4);
# 
#   for (t in 2:T)
#     y[t] ~ normal(y[t-1] + beta[period[t]], sigma_proc);
# 
#   for (i in 1:N_obs)
#     y_obs[i] ~ normal(y[obs_t[i]], sigma_obs);
# }
# "

# stan_code <- "
# data {
#   int<lower=1> T;          // number of years
#   int<lower=1> P;          // number of periods
#   int<lower=1> N_obs;      // number of observations
# 
#   int<lower=1,upper=T> obs_t[N_obs];  // indices of observations
#   vector[N_obs] y_obs;                // observed values
# 
#   int<lower=1,upper=P> period[T];     // period assignment for each year
# 
#   vector[P] beta_mu;                  // prior mean for beta
#   vector[P] beta_se;                  // prior sd for beta
# }
# parameters {
#   vector[T] eta;        // standard normal latent increments
#   vector[P] beta;       // period effects
# 
#   real<lower=0> sigma_proc;  // process noise
#   real<lower=0> sigma_obs;   // observation noise
# }
# transformed parameters {
#   vector[T] y;         // latent state
# 
#   // initialize first year
#   y[1] = 355 + 8.4 * eta[1];
# 
#   // random-walk for subsequent years
#   for (t in 2:T)
#     y[t] = y[t-1] + beta[period[t]] + sigma_proc * eta[t];
# }
# model {
#   // standard normal prior for eta (non-centered)
#   eta ~ normal(0, 1);
# 
#   // priors
#   beta ~ normal(beta_mu, beta_se);
#   //sigma_proc ~ normal(0, 1);
#   //sigma_obs ~ normal(0, 1);
#   //sigma_proc ~ normal(0.5, 0.1) T[0, ];
#   sigma_obs ~ exponential(1);
# 
#   // likelihood
#   y_obs ~ normal(y[obs_t], sigma_obs);
# }
# "

stan_code <- "
data {
  int<lower=1> T;
  int<lower=1> P;
  int<lower=1> N_obs;

  int<lower=1,upper=T> obs_t[N_obs];
  vector[N_obs] y_obs;

  int<lower=1,upper=P> period[T];

  vector[P] beta_mu;
  vector[P] beta_se;
}
parameters {
  real y0;                 // initial absolute level
  vector[P] beta;          
  real<lower=0> sigma_obs;
}
transformed parameters {
  vector[T] y;

  y[1] = y0;

  for (t in 2:T)
    y[t] = y[t-1] + beta[period[t]];
}
model {
  // anchor absolute level with prior
  y0 ~ normal(355, 8.4);

  beta ~ normal(beta_mu, beta_se);
  sigma_obs ~ exponential(1);

  y_obs ~ normal(y[obs_t], sigma_obs);
}
"


years <- 1990:2025
T <- length(years)

period <- cut(
  years,
  breaks = c(1990, 2000, 2010, 2020, 2025),
  labels = 1:4,
  # labels = paste0('p', 1:4),
  right = TRUE, 
  include.lowest = TRUE
) |> as.integer()

##---------------------below are the things I have to update--------------------------

beta_mu <- c(1.81, 2.04, 2.36, 2.85)
beta_se <- c(0.5, 0.5, 0.6, 0.6)


# obs <- data.frame(
#   year  = c(1993,2000, 2007, 2012, 2018),
#   value = c(355, 368, 380, 397, 410)
# )

obs <- data.frame(
    year  = c(1993, 2007, 2018),
    value = c(355, 380, 410)
  )

obs_t <- match(obs$year, years)

stan_data <- list(
  T = T,
  P = 4,
  N_obs = nrow(obs),
  obs_t = obs_t,
  y_obs = obs$value,
  period = period,
  beta_mu = beta_mu,
  beta_se = beta_se
)


library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

fit <- stan(
  model_code = stan_code,
  data = stan_data,
  chains = 4,
  iter = 2000,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 20
  )
)

summary(fit)

# get y values and beta
post <- rstan::extract(fit)

y_mean <- apply(post$y, 2, mean)
y_ci   <- apply(post$y, 2, quantile, probs = c(0.025, 0.975))

df <- data.frame(
  year = years,
  mean = y_mean,
  lo = y_ci[1,],
  hi = y_ci[2,]
)






