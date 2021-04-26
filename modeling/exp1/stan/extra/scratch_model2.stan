// Title: Univariate Normal Non-linear Exepected Utility Maximization Model
// Author: Tyler Adkins
// Date: June 22nd, 2020

functions {
  int get_max_idx(real[] EU);
  #include "/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/stan/hpp/get_max_idx.hpp"
  real get_aim(vector grid,  // grid of aim-points
               vector p_gain,  // probability of gain, given aim-point
               vector p_loss,  // probability of loss, given aim-point
               real gain,  // gain value
               real loss,  // loss value
               real beta,  // loss decay
               real lambda  // loss aversion
               ) {
    real EU[rows(grid)]; // container for expected utilities
    // compute expected utility of each aim point
    EU = to_array_1d(p_gain * gain + p_loss * (-lambda*(-loss)^beta));
    return grid[get_max_idx(EU)];  // return the aimpoint with max EU
  }
}

data {
  int<lower=0> N;  // number of observations/trials
  int<lower=1> G;  // number of aimpoints to consider
  int<lower=1> C;  // number of conditions
  int<lower=1> S;  // number of individuals
  int<lower=1> sub_id[N];  // individiual of each trial
  int<lower=1> cond_id[N];  // condition of each trial
  vector[N] gain;  // gain of each trial
  vector[N] loss;  // loss of each trial
  vector[N] obs;  // observed x-endpoint of each trial  
  vector[C] cond;  // list of unique conditions (combinations of gain & loss)
  vector[G] grid;  // candidate aim-points
  matrix<lower=0,upper=1>[G,S] p_gain;  // gain probabilities for each aim-point
  matrix<lower=0,upper=1>[G,S] p_loss;  // loss probabilities for each aim-point
  int<lower=0,upper=1> get_obs_rep;  // get posterior predicted observations?
  int<lower=0,upper=1> get_log_lik;  // get log-liklihoods?
  int<lower=0,upper=1> prior_only;  // ignore the observations?
}
parameters {
  real beta_mu;
  real<lower=0> beta_sd;
  vector[S] beta_re;
  // real lambda_mu;
  // real<lower=0> lambda_sd; 
  // vector[S] lambda_re;
  real<lower=0> sigma_mu;
  real<lower=0> sigma_sd;
  vector[S] sigma_re;
}
transformed parameters {
  vector[S] beta;  // loss decay
  vector[S] lambda;  // loss aversion
  vector<lower=0>[S] sigma;  // individual error
  matrix[C,S] mu_cond;  // condition/subject specific aimpoints
  
  beta = beta_mu + beta_sd*beta_re;
  lambda = lambda_mu + lambda_sd*beta_re;
  sigma = sigma_mu + sigma_sd*sigma_re;
  
  // pre-compute condition/subject means
  for (c in 1:C) {
    for (s in 1:S) {
      mu_cond[c,s] = lambda[s] + cond[c]*beta[s];
      // get_aim(grid, p_gain[,s], p_loss[,s], 1, cond[c], beta[s], lambda[s]);
    }
  }
}
model {
  beta_mu ~ normal(1,.5);
  beta_sd ~ normal(.5,.25);
  beta_re ~ std_normal();
  lambda_mu ~ normal(1,.5);
  lambda_sd ~ normal(.5,.25);
  lambda_re ~ std_normal();
  sigma_mu ~ normal(.5,.25);
  sigma_sd ~ normal(.5,.25);
  sigma_re ~ std_normal();

  if (prior_only == 0) {
    {
      vector[N] Mu;
      vector[N] Sigma;
      // organize means
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
