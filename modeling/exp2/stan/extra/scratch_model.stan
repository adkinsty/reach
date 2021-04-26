// Title: Univariate Normal Non-linear Exepected Utility Maximization Model
// Author: Tyler Adkins
// Date: June 22nd, 2020

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
  int<lower=1> C;  // number of conditions
  int<lower=1> S;  // number of individuals
  int<lower=2,upper=2> K; // number of utility function parameters
  int<lower=1> sub_id[N];  // individiual of each trial
  int<lower=1> cond_id[N];  // condition of each trial
  vector[C] cond; // list of unique conditions (combinations of gain & loss)
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
  cholesky_factor_corr[K] Lcorr; // prior cholesky factor
  vector<lower=0>[K] tau; // prior scale
  vector<lower=0>[K] theta_mu;  // group mean coeffs
  vector<lower=0>[K] theta[S]; // indivual coeffs
  real<lower=0> sigma_mu;  // group mean error
  real<lower=0> sigma_sd;  // sd of error across individuals
  vector<lower=0>[S] sigma;  // individual error
}
transformed parameters {
  matrix[C,S] mu_cond;  // condition/subject specific aimpoints
  // pre-compute condition/subject means
  for (c in 1:C) {
    for (s in 1:S) {
      mu_cond[c,s] = get_aim(grid, p_gain[,s], p_loss[,s], 1, cond[c], theta[s][1], theta[s][2]);
    }
  }
}
model {
  tau ~ normal(0, 1);
  Lcorr ~ lkj_corr_cholesky(2);
  to_vector(theta_mu) ~ normal(1,.5);
  theta ~ multi_normal_cholesky(theta_mu,diag_pre_multiply(tau,Lcorr));
  
  sigma_mu ~ normal(0.5,0.1);
  sigma_sd ~ normal(0.05,0.01);
  sigma ~ normal(sigma_mu,sigma_sd);
  
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
