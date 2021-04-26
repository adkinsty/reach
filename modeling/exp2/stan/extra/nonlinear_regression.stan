// Title: Hierarchical linear regression model
// Description: 
//   -optimized by scaling, re-parameterizing, and vectorizing
//   -regression coefficients and sigma vary by subject
// Author: Tyler Adkins
// Date: June 24th, 2020

data {
  int<lower=0> N;                    // num trials (total)
  int<lower=2,upper=2> K;            // num trial predictors
  int<lower=1> J;                    // num subjects
  int<lower=1> L;                    // num subject predictors
  int<lower=1,upper=3> C;            // num conditions
  vector<lower=-5,upper=0>[C] cond;  // unique conditions
  int<lower=1> cc[N];                // cond idx for trial
  int<lower=1> jj[N];                // subject idx for trial
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
  matrix[K, J] re;                            // random effect on beta
  vector<lower=0>[K] tau;                     // prior scale on re
  cholesky_factor_corr[K] L_Omega;            // factor for re corr
  matrix[L, K] gamma;                // group coeffs (for betas)
  real<lower=0> sigma_m;                      // group mean sigma
  real<lower=0> sigma_s;                      // sd of sigma between subjects
  vector[J] sigma_re;                         // individual random effect on sigma
}
transformed parameters {
  matrix[J, K] beta;
  vector<lower=0>[J] sigma;
  matrix[J,C] mu_cond;
  beta = u * gamma + (diag_pre_multiply(tau,L_Omega) * re)'; 
  sigma = sigma_m + sigma_s * sigma_re;

  for (j in 1:J) {
    for (c in 1:C) {
      mu_cond[j,c] =  beta[j,1] + beta[j,2]*cond[c];    
    }
  }
}
model {
  L_Omega ~ lkj_corr_cholesky(2);
  tau ~ std_normal();
  to_vector(re) ~ std_normal();
  to_vector(gamma) ~ beta(//std_normal();
  
  sigma_m ~ std_normal();
  sigma_s ~ std_normal();
  sigma_re ~ std_normal();
  
  {
    vector[N] Mu;
    for (n in 1:N) {
      Mu[n] = mu_cond[jj[n],cc[n]];
    }
   zy ~ normal(Mu, sigma[jj]);
  }
}
generated quantities {
  matrix[K,K] Omega;
  vector[N] zyrep;
  vector[N] yrep;
  real log_lik[N];
  
  Omega = multiply_lower_tri_self_transpose(L_Omega);

  if (get_yrep == 1) {
    for (n in 1:N) {
      zyrep[n] = normal_rng(dot_product(beta[jj[n]] , x[n,]), sigma[jj[n]]);
      yrep = zyrep * sample_sd + sample_mean;
    }
  }
  if (get_log_lik == 1) {
    for (n in 1:N) {
      log_lik[n] = normal_lpdf(zy[n] | dot_product(beta[jj[n]] , x[n,]), 
                                       sigma[jj[n]]);
    }
  }
}
