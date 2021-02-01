functions {
  real switch_eta(real t, real t1, real eta, real nu, real xi) {
    return(eta + (1 - eta) / (1 + exp(xi * (t - t1 - nu))));
  }
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      int N = x_i[1];
      real tswitch = x_r[1];
      
      real beta = theta[1];
      real gamma = theta[2];
      real a = theta[3];
      real eta = theta[4]; // reduction in transmission rate after quarantine
      real nu = theta[5]; // shift of quarantine implementation
      real xi = theta[6]; // slope of quarantine implementation
      real i0 = theta[7]; // initial number of infected people
      real e0 = theta[8]; // initial number of infected people
      real forcing_function = switch_eta(t,tswitch,eta,nu,xi); // switch function
      real beta_eff = beta * forcing_function; // beta decreased to take control measures into account
      real init[4] = {N - i0 - e0, e0, i0, 0}; // initial values

      real S = y[1] + init[1];
      real E = y[2] + init[2];
      real I = y[3] + init[3];
      real R = y[4] + init[4];
      
      real dS_dt = -beta_eff * I * S / N;
      real dE_dt =  beta_eff * I * S / N - a * E;
      real dI_dt = a * E - gamma * I;
      real dR_dt =  gamma * I;
      
      return {dS_dt, dE_dt, dI_dt, dR_dt};
  }
}
data {
  int<lower=1> n_days;
  real t0;
  real tswitch; // date of introduction of control measures
  real ts[n_days];
  int N; // population size
  int cases[n_days];
  int t_survey_start; // antibody survey data
  int t_survey_end;
  int n_infected_survey;
  int n_tested_survey;
}
transformed data {
  int x_i[1] = { N }; //formatting to feed the ODE function
  real x_r[1] = {tswitch};
}
parameters {
  real<lower=0> gamma; // SEIR parameters
  real<lower=0> beta;
  real<lower=0> a;
  real<lower=0> phi_inv; 
  real<lower=0,upper=1> eta; // reduction in transmission due to control measures (in proportion of beta)
  real<lower=0> nu; // shift of quarantine implementation (strictly positive as it can only occur after tswitch)
  real<lower=0,upper=1> xi_raw; // slope of quarantine implementation (strictly positive as the logistic must be downward)
  real<lower=0> i0; // number of infected people inititally
  real<lower=0> e0; // number of exposed people inititally
  real<lower=0, upper=1> p_reported; // probability for an infected person to be reported (i.e counted as a case)
}
transformed parameters{
  real y[n_days, 4];
  real incidence[n_days - 1];
  real phi = 1. / phi_inv;
  real xi = xi_raw + 0.5;
  real theta[8];
  real<lower=0, upper=1> p_infected_survey; //proportion of people having been infected at week 5 (between 4 and 7 may)
  theta = {beta, gamma, a, eta, nu, xi, i0, e0};
  y = integrate_ode_rk45(sir, rep_array(0.0, 4), t0, ts, theta, x_r, x_i);
  for (i in 1:n_days-1){
    incidence[i] = -(y[i+1, 2] - y[i, 2] + y[i+1, 1] - y[i, 1]) * p_reported; //-(E(t+1) - E(t) + S(t+1) - S(t))
  }
  // mean number of infected + recovered people during week 5
  p_infected_survey = mean(to_vector(y[t_survey_start:t_survey_end, 4])) / N;
}
model {
  //priors
  beta ~ normal(2, 1);
  gamma ~ normal(0.4, 0.5);
  a ~ normal(0.4, 0.5);
  phi_inv ~ exponential(5);
  i0 ~ normal(0, 10);
  e0 ~ normal(0, 10);
  p_reported ~ beta(1, 2);
  eta ~ beta(2.5, 4);
  nu ~ exponential(1./5);
  xi_raw ~ beta(1, 1);
  
  //sampling distribution
  n_infected_survey ~ binomial(n_tested_survey, p_infected_survey); // we fit the survey data to our latent parameter
  //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people 
  cases[1:(n_days-1)] ~ neg_binomial_2(incidence, phi);
}
generated quantities {
  real R0 = beta / gamma;
  real Reff[n_days]; // R0 but taking into account environmental changes
  real recovery_time = 1 / gamma;
  real incubation_time = 1 / a;
  real pred_cases[n_days-1];
  pred_cases = neg_binomial_2_rng(incidence, phi);
  for (i in 1:n_days)
    Reff[i] = switch_eta(i, tswitch, eta, nu, xi) * beta / gamma;
}

