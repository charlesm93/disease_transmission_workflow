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

generated quantities {
    real gamma;
    real beta;
    real phi;
    real phi_inv;
    real R0;
    real recovery_time;
    real theta[2];
    real pred_cases[T];
    real y[T,3];
    phi_inv = fabs(cauchy_rng(0.5, 1));
    phi = 1. / phi_inv;
    gamma = exponential_rng(3.);
    beta = fabs(normal_rng(0, 0.4));
    R0 = beta/gamma;
    recovery_time = 1. / gamma;
    theta[1] = beta;
    theta[2] = gamma;
    y = integrate_ode_bdf(sir, y0, t0, ts, theta, x_r, x_i);
    pred_cases = neg_binomial_2_rng(col(to_matrix(y),2) + 1e-9, phi);
}
