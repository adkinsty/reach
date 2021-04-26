// Title: Linear regression model
// Author: Tyler Adkins
// Date: June 24th, 2020

data {
  int<lower=0> N;                    // num trials (total)
  int<lower=1> K;                    // num trial predictors
  int<lower=1> J;                    // num subjects
  int<lower=1> L;                    // num subject predictors
  int<lower=1> jj[N];                // subject for trial
  int<lower=1> C;                    // number of conditions
  vector[C] cond;                    // conditions
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
  vector[N] zy = (y - sample_mean)/sample_sd;   // z-score outcomes
}
parameters {
  matrix[K, J] beta_re;            
  vector<lower=0>[K] beta_tau;       
  cholesky_factor_corr[K] beta_LO;   
  matrix[L, K] beta_gamma;           
  real<lower=0> sigma_m;  
  real<lower=0> sigma_s;
  vector<lower=0>[J] sigma;
  real theta; 
}
transformed parameters {
  matrix[J, K] beta;
  vector<lower=0>[N] sigma_c; 

  beta = u * beta_gamma + (diag_pre_multiply(beta_tau,beta_LO) * beta_re)'; 
  
  for (n in 1:N) {
    sigma_c[n] = sigma[jj[n]] + theta*x[n,2];
  }
  // non-centered parameterizations for:
  //   beta ~ multi_normal_cholesky(beta_mu, beta_corr)
}
model {
  beta_LO ~ lkj_corr_cholesky(2);
  beta_tau ~ std_normal();
  to_vector(beta_re) ~ std_normal();
  to_vector(beta_gamma) ~ std_normal();
  
  sigma_m ~ std_normal();
  sigma_s ~ std_normal();
  sigma ~ std_normal();
  theta ~ std_normal();
  
  if (prior_only == 0) {
    zy ~ normal(rows_dot_product(beta[jj] , x), sigma);
  }
}
generated quantities {
  vector[N] zyrep;
  vector[N] yrep;
  real log_lik[N];

  if (get_yrep == 1) {
    for (n in 1:N) {
      zyrep[n] = normal_rng(dot_product(beta[jj[n]] , x[n,]), dot_product(sigma[jj[n]] , x[n,]));
      yrep[n] = zyrep[n] * sample_sd + sample_mean;
    }
  }
  if (get_log_lik == 1) {
    for (n in 1:N) {
      log_lik[n] = normal_lpdf(zy[n] | dot_product(beta[jj[n]], x[n,]), dot_product(sigma[jj[n]] , x[n,]));
    }
  }
}
