// Title: Univariate Normal Non-linear Exepected Utility Maximization Model
// Author: Tyler Adkins
// Date: June 22nd, 2020

data {
  int<lower=0> N;                    // num trials (total)
  int<lower=1> K;                    // num trial predictors
  int<lower=1> J;                    // num subjects
  int<lower=1> L;                    // num subject predictors
  int<lower=1> jj[N];                // subject for trial
  matrix[N, K] x;                    // trial predictors
  matrix[J, L] u;                    // subject predictors
  vector[N] y;                       // outcomes
  int<lower=0,upper=1> get_yrep;     // get posterior predicted observations?
  int<lower=0,upper=1> get_log_lik;  // get log-liklihoods?
  int<lower=0,upper=1> prior_only;   // ignore the observations?
}
transformed data {
  real sample_mean = mean(y);
  real sample_sd = sd(y);
  vector[N] zy = (y - sample_mean)/sample_sd;   // z-score
}
parameters {
  matrix[K, J] re;                             // random effect on beta
  vector<lower=0>[K] tau;                     // prior scale on re
  cholesky_factor_corr[K] L_Omega;            // factor for re corr
  matrix[L, K] gamma;                         // group coef
  real<lower=0> sigma;
}
transformed parameters {
  matrix[J, K] beta;
  beta = u * gamma + (diag_pre_multiply(tau,L_Omega) * re)'; 
  // non-centered parameterizations for:
  // beta:
  //   beta ~ multi_normal_cholesky(beta_mu, beta_corr)
  // sigma:
  //   sigma ~ normal(sigma_m, sigma_s)
}
model {
  L_Omega ~ lkj_corr_cholesky(2);
  tau ~ std_normal();
  to_vector(re) ~ std_normal();
  to_vector(gamma) ~ std_normal();
  
  sigma ~ std_normal();

  zy ~ normal(rows_dot_product(beta[jj] , x), sigma);
}
generated quantities {
  matrix[K,K] Omega;
  vector[N] zyrep;
  vector[N] yrep;
  real log_lik[N];
  matrix[L,K] beta_mu;

  beta_mu = u * gamma;

  Omega = multiply_lower_tri_self_transpose(L_Omega);

  if (get_yrep == 1) {
    for (n in 1:N) {
      zyrep[n] = normal_rng(dot_product(beta[jj[n]] , x[n,]), sigma);
      yrep = zyrep * sample_sd + sample_mean;
    }
  }
  if (get_log_lik == 1) {
    for (n in 1:N) {
      log_lik[n] = normal_lpdf(zy[n] | dot_product(beta[jj[n]] , x[n,]), sigma);
    }
  }
}
