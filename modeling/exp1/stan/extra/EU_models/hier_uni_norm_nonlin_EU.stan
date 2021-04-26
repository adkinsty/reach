// Title: Univariate Normal Non-linear Exepected Utility Maximization Model
// Author: Tyler Adkins
// Date: June 25th, 2020

functions {
  int get_max_idx(real[] EU);
  #include "/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/stan/hpp/get_max_idx.hpp"
  real get_aim(vector grid, vector p_gain, vector p_loss, real gain, 
                 real loss, real beta, real lambda) {
    // Get aimpoint that maximizes expected utility
    real EU[rows(grid)];
    EU = to_array_1d(p_gain * gain + p_loss * (-lambda*(-loss)^beta));
    return grid[get_max_idx(EU)];
  }
}

data {
  int<lower=0> N; // number of observations
  int<lower=1> G; // number of aimpoints to consider
  int<lower=1> C;  // number of unique conditions
  int<lower=1> S;  // number of unique subjects
  int<lower=1> sub_id[N];  // subject IDs
  int<lower=1> cond_id[N];  // condition IDs (1,2,3)
  vector[C] cond; // list of unique conditions
  vector[N] gain; // gain covariate
  vector[N] loss; // gain covariate
  vector[G] grid; // candidate aim-points
  matrix<lower=0,upper=1>[G,S] p_gain; // gain probabilities for each aim-point
  matrix<lower=0,upper=1>[G,S] p_loss; // loss probabilities for each aim-point
  vector[N] obs; // observations  
  // generated quantities, y/n
  int<lower=0,upper=1> get_obs_rep;
  int<lower=0,upper=1> get_log_lik;
  int<lower=0,upper=1> prior_only;  
}
parameters {
  real<lower=0> beta_m;
  real<lower=0> beta_s;
  vector<lower=0>[S] beta_re;
  real<lower=0> lambda_m;
  real<lower=0> lambda_s;
  vector<lower=0>[S] lambda_re;
  real<lower=0> sigma_m;  
  real<lower=0> sigma_s; 
  vector<lower=0>[S] sigma_re;
  vector<lower=0>[S] bias;  // individual bias
}
transformed parameters {
  matrix[C,S] mu_cond;  // condition/subject specific aimpoints
  vector<lower=0>[S] beta;  
  vector<lower=0>[S] lambda;  
  vector<lower=0>[S] sigma;
  beta = to_vector(beta_m + beta_re*beta_s);
  lambda = lambda_m + lambda_re*lambda_s;
  sigma = sigma_m + sigma_re*sigma_s;
  // pre-compute condition/subject means
  for (c in 1:C) {
    for (s in 1:S) {
      mu_cond[c,s] = bias[s] + get_aim(grid, p_gain[,s], p_loss[,s], 1, cond[c], beta[s], lambda[s]);
    }
  }
}
model {
  // group
  beta_m ~ gamma(2,2);
  beta_s ~ gamma(2,2);

  lambda_m ~ gamma(4,2);
  sigma_m ~ gamma(2,3);
  
  bias ~ normal(0,.5);
  
  if (prior_only == 0) {
    {
      vector[N] Mu;
      vector[N] Sigma;
      // set up predictors
      for (n in 1:N) {
        Mu[n] = mu_cond[cond_id[n],sub_id[n]];
        Sigma[n] = sigma[sub_id[n]];
      }
      // get likelihood
      obs ~ normal(Mu,Sigma);
    }
  }
}
generated quantities {
    vector[N] obs_rep;
    vector[N] log_lik;
    
    if (get_obs_rep == 1) {
        for (n in 1:N) {
            obs_rep[n] = normal_rng(mu_cond[cond_id[n],sub_id[n]],sigma[sub_id[n]]);
        }    
    }
    if (get_log_lik == 1) {
        for (n in 1:N) {
            log_lik[n] = normal_lpdf(obs[n] | mu_cond[cond_id[n],sub_id[n]],sigma[sub_id[n]]);
        }
    }
}
