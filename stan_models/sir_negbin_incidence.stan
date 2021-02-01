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
    real dydt[3];
    dydt[1] =  -  theta[1] * y[1] * y[2] / x_i[1];
    dydt[2] =    theta[1] * y[1] * y[2] / x_i[1] - theta[2] * y[2];
    dydt[3] = theta[2] * y[2];
    return dydt;
  }
}

data {
  int<lower=1> T;
  real y0[3];
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
  real<lower=0> phi_inv;
}

transformed parameters{
  real y[T,3];
  real incidence[T];
  real phi = 1. / phi_inv;
  {
    real theta[2];
    theta[1] = beta;
    theta[2] = gamma;

    y = integrate_ode_bdf(sir, y0, t0, ts, theta, x_r, x_i);
    incidence[1] = y[1, 2];
    for (i in 2:T){
      incidence[i] = y[i, 2] - y[i-1, 2] + y[i, 3] - y[i-1, 3];
      }
  }
}


model {
  beta ~ lognormal(0, 1);
  gamma ~ lognormal(0, 1);
  
  //col(matrix x, int n) - The n-th column of matrix x
  phi_inv ~ cauchy(0.5, 1);
  cases ~ neg_binomial_2(incidence, phi);

}

generated quantities {
  real R0 = beta/gamma;
  real pred_cases[T];
  pred_cases = neg_binomial_2_rng(incidence, phi);
}
