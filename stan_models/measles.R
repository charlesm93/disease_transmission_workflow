library(outbreaks)
library(rstan)
library(dplyr)
rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())
set.seed(5) # for reproductibility



df_confirmed <- read.csv("data/niamey.csv")
df_confirme
View(df_confirmed)

df_confirmed %>% filter(community=="B") %>% ggplot() + geom_line(mapping=aes(x=biweek, y=measles))

# time series of cases
cases <- (df_confirmed %>% filter(community=="B"))$measles  #Number of students in bed

# total count
N <- 1000;

# times
T <- length(cases); t <- seq(0, T, by = 1); t0 = 0; t <- t[-1]; 

#initial conditions
i0 <- 1; s0 <-  N-i0; r0 <- 0; y0 = c( S=s0, I=i0, R=r0);

# data for Stan
data_sir <-  list (T  = T, y0 = y0, t0 = t0, ts = t, N=N, cases=cases)

# number of MCMC steps
niter <- 1000

?influenza_england_1978_school
model <- stan_model("sir_negbin_incidence.stan")

fit_sir_negbin <- sampling(model, 
                           data = data_sir, 
                           iter=niter,
                           chains=2,
                           refresh=20,
                           control=list(adapt_delta=0.99))


pars=c('beta', 'gamma', 'phi_inv')

summary(fit_sir_negbin, pars=pars, probs=c(0.025, 0.5, 0.975))$summary

stan_dens(fit_sir_negbin,pars=pars)

stan_dens(fit_sir_negbin,pars=pars, separate_chains = TRUE)


s <- rstan::extract(fit_sir_negbin)

I_est <- s$pred_cases
df_est = data.frame(  I_median=apply(I_est, 2, median), 
                      I_low=apply(I_est, 2, quantile, probs=c(0.025)),
                      I_high=apply(I_est, 2, quantile, probs=c(0.975)), 
                      t)

ggplot(df_est, mapping=aes(x=t))+
  geom_point(mapping=aes(y=cases)) +
  geom_ribbon(aes(ymin = I_low, ymax = I_high), fill = "orange", alpha = 0.6) + 
  geom_line(data=df_est, mapping=aes(x=t, y=I_median)) +
  labs(x = "Day", y="Number of students in bed") +
  labs(colour = "Median+95%CI")

