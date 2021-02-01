functions {
  //theta[1] = beta
  //theta[2] = gamma
  //x_i[1] = N, total population
  // y = S, I, R
  real[] sir(real t,
             real[] y,
             real[] theta,
             real[] x_r,
             int[] x_i) {
    real dydt[4];
    real beta = theta[1];
    real gamma = theta[2];
    real a = theta[3];
    real S = y[1];
    real E = y[2];
    real I = y[3];
    real R = y[4];
    int N = x_i[1];
    dydt[1] =  -  beta * S * I / N;
    dydt[2] =    beta * S * I / N - a * E;
    dydt[3] = a * E - gamma * I;
    dydt[4] = gamma * I; 
    return dydt;
  }
}

data {
  int<lower=1> T;
  real y0[4];
  real t0;
  real ts[T];
  int N;
  int cases[T];
}

transformed data {
  real x_r[0];
  int x_i[1];
  x_i[1]=N;
}

parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> a;
  real phi_inv;
}

transformed parameters{
  real y[T,4];
  real incidence[T];
  real phi = 1. / phi_inv;
  {
    real theta[3];
    theta[1] = beta;
    theta[2] = gamma;
    theta[3] = a;

    y = integrate_ode_bdf(sir, y0, t0, ts, theta, x_r, x_i);
    incidence[1] = y[1, 3];
    for (i in 2:T){
      incidence[i] = y[i, 3] - y[i-1, 3] + y[i, 4] - y[i-1, 4];
      }
  }
}


model {
  beta ~ lognormal(0, 1);
  gamma ~ lognormal(0, 1);
  a ~ lognormal(0, 1);
  //col(matrix x, int n) - The n-th column of matrix x
  phi_inv ~ cauchy(0.5, 1);
  cases ~ neg_binomial_2(incidence, phi);

}

generated quantities {
  real R0 = beta/gamma;
  real pred_cases[T];
  pred_cases = neg_binomial_2_rng(incidence, phi);
}
