functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real N = x_i[1];
      
      real beta = theta[1];
      real gamma = theta[2];
      real a = theta[3];
      real i0 = theta[4];
      real e0 = theta[5];
      
      real init[4] = {N - i0 - e0, e0, i0, 0}; // initial values
      real S = y[1] + init[1];
      real E = y[2] + init[2];
      real I = y[3] + init[3];
      real R = y[4] + init[4];
      
      real dS_dt = -beta * I * S / N;
      real dE_dt =  beta * I * S / N - a * E;
      real dI_dt = a * E - gamma * I;
      real dR_dt =  gamma * I;
      
      return {dS_dt, dE_dt, dI_dt, dR_dt};
  }
}
data {
  int<lower=1> n_days;
  real t0;
  real ts[n_days];
  int N;
  int cases[n_days];
}
transformed data {
  real x_r[0];
  int x_i[1] = { N };
}
parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> inv_a;
  //real<lower=0> inv_a;
  real<lower=0> phi_inv;
  real<lower=0, upper=1> p_reported; // proportion of infected (symptomatic) people reported
  real<lower=0> i0;
  real<lower=0> e0;
}
transformed parameters{
  real y[n_days, 4];
  real incidence[n_days - 1];
  real phi = 1. / phi_inv;
  real a = 1. / inv_a;
  real theta[5] = {beta, gamma, a, i0, e0};

  //y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  y = integrate_ode_rk45(sir, rep_array(0.0, 4), t0, ts, theta, x_r, x_i);
  
  for (i in 1:n_days-1){
    incidence[i] = -(y[i+1, 2] - y[i, 2] + y[i+1, 1] - y[i, 1]) * p_reported; //-(E(t+1) - E(t) + S(t+1) - S(t))
  }
}
model {
  //priors
  beta ~ normal(2, 1);
  gamma ~ normal(0.4, 0.5);
  //a ~ normal(0.1, 0.05);
  inv_a ~ normal(6, 1);
  //a ~ gamma(2, 1./20);
  //inv_a ~ normal(10, 2);//normal(7, 3);
  phi_inv ~ exponential(5);
  p_reported ~ beta(1, 2);
  i0 ~ normal(0, 10);
  e0 ~ normal(0, 10);

  //sampling distribution
  //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people 
  cases[1:(n_days-1)] ~ neg_binomial_2(incidence, phi);
}
generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  real incubation_time = 1 / a;
  real pred_cases[n_days-1];
  pred_cases = neg_binomial_2_rng(incidence, phi);
}

