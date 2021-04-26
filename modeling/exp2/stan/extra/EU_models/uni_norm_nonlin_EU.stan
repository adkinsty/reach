// Title: Univariate Normal Nonlinear Expected Utility Model
// Author: Tyler Adkins
// Date: June 19th, 2020

functions {
  int get_max_idx(real[ ] EU);
  #include "/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/stan/hpp/get_max_idx.hpp"
  real u(real v, real beta, real lambda){
    // Utility function
    if (v >= 0) {
      return v;
    } else {
      return -lambda * (-v)^beta;
    }
  }
  real get_aim(vector grid, vector p_gain, vector p_loss, real gain, 
               real loss, real beta, real lambda) {
    // Get aimpoint that maximizes expected utility
    real aim;
    real EU[rows(grid)];
    EU = to_array_1d(p_gain * u(gain,beta,lambda) + p_loss * u(loss,beta,lambda));
    aim = grid[get_max_idx(EU)];
    return aim;
  }
}

data {
  int<lower=0> N; // number of observations
  int<lower=1> G; // number of aimpoints to consider
  int<lower=1> C; // number of unique conditions
  int<lower=1> cond_id[N];  // condition IDs (1,2,3)
  vector[C] cond; // list of unique conditions
  vector[G] grid; // candidate aim-points
  vector[G] p_gain; // gain probabilities for each aim-point
  vector[G] p_loss; // loss probabilities for each aim-point
  vector[N] gain; // gain covariate
  vector[N] loss; // gain covariate
  vector[N] obs; // observations  // generated quantities, y/n
  int<lower=0> get_obs_rep;
  int<lower=0,upper=1> get_log_lik;
  int<lower=0,upper=1> prior_only;  
}
parameters {
  real<lower=0> beta;   // loss decay parameter
  real<lower=0> lambda; // loss aversion parameter
  real<lower=0> sigma;  // endpoint variance
}
transformed parameters {
  vector[C] mu_cond;
  for (c in 1:C) {
    mu_cond[c] = get_aim(grid, p_gain, p_loss, 1, cond[c], beta, lambda);
  }
}
model {
  // priors
  beta ~ normal(1,.5);
  lambda ~ normal(1,.5);
  sigma ~ normal(.5,.5);
  // likelihood
  if (prior_only == 0) {
    vector[N] mu;
    for (n in 1:N) {
      mu[n] = mu_cond[cond_id[n]];
    }
    obs ~ normal(mu, sigma);
  }
}
generated quantities {
    vector[N] obs_rep;
    vector[N] log_lik;
    
    if (get_obs_rep == 1) {
        for (n in 1:N) {
            obs_rep[n] = normal_rng(mu_cond[cond_id[n]],sigma);
        }    
    }
    if (get_log_lik == 1) {
        for (n in 1:N) {
            log_lik[n] = normal_lpdf(obs[n] | mu_cond[cond_id[n]],sigma);
        }
    }
}
