// Title: Multivariate Normal Non-linear Exepected Utility Maximization Model
// Author: Tyler Adkins
// Date: June 22nd, 2020

functions {
  int get_max_idx(real[] EU);
  #include "/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/stan/hpp/get_max_idx.hpp"
  real u(real v, real beta, real lambda){
    // Utility function
    if (v >= 0) {
      return v;
    } else {
      return -lambda * (-v)^beta;
    }
  }
  vector get_aim(vector grid, real[] p_gain, real[] p_loss, real gain, 
                 real loss, real beta, real lambda) {
    // Get aimpoint that maximizes expected utility
    vector[2] aim;
    real EU[rows(grid)];
    EU = to_array_1d(to_vector(p_gain) * u(gain,beta,lambda) + to_vector(p_loss) * u(loss,beta,lambda));
    aim[1] = grid[get_max_idx(EU)]; // maximizing x aimpoint
    aim[2] = 0; // never aim above midline
    return aim;
  }
}

data {
  int<lower=0> N; // number of observations
  int<lower=1> G; // number of aimpoints to consider
  int<lower=2,upper=2> J;  // number of outcome-dim
  int<lower=1> C;  // number of unique conditions
  int<lower=1> S;  // number of unique subjects
  int<lower=1> sub_id[N];  // subject IDs
  int<lower=1> cond_id[N];  // condition IDs (1,2,3)
  vector[C] cond; // list of unique conditions
  vector[N] gain; // gain covariate
  vector[N] loss; // gain covariate
  vector[G] grid; // candidate aim-points
  real<lower=0,upper=1> p_gain[G,S]; // gain probabilities for each aim-point
  real<lower=0,upper=1> p_loss[G,S]; // loss probabilities for each aim-point
  vector[J] obs[N]; // observations  
  // generated quantities, y/n
  int<lower=0,upper=1> get_obs_rep;
  int<lower=0,upper=1> get_log_lik;
  int<lower=0,upper=1> prior_only;  
}
parameters {
  real<lower=0> m_beta;  // group mean beta
  real<lower=0> m_lambda;  // group mean lambda
  real<lower=0> s_beta;  // group sd beta
  real<lower=0> s_lambda;  // group sd lambda
  real<lower=0> beta[S];   // loss decay parameter
  real<lower=0> lambda[S]; // loss aversion parameter
  vector<lower=0>[J] sigma;
  cholesky_factor_corr[J] L_corr;
}
transformed parameters {
  vector[J] mu_cond[C,S];
  for (c in 1:C) {
    for (s in 1:S) {
      mu_cond[c,s] = get_aim(grid, p_gain[,s], p_loss[,s], 1, cond[c], beta[s], lambda[s]);
    }
  }
}
model {
  m_beta ~ normal(1,0.25);
  m_lambda ~ normal(1,.25);
  s_beta ~ normal(.5,.25);
  s_lambda ~ normal(.5,.25);
  beta ~ normal(m_beta,s_beta);
  lambda ~ normal(m_lambda,s_lambda);
  sigma ~ normal(0,1);
  L_corr ~ lkj_corr_cholesky(2.0);
  // likelihood
  if (prior_only == 0) {
    {
      vector[J] mu[N];
      for (n in 1:N) {
        mu[n] = mu_cond[cond_id[n],sub_id[n]];
      }
      obs ~ multi_normal_cholesky(mu, diag_pre_multiply(sigma, L_corr));
    }
  }
}
generated quantities {
    corr_matrix[J] corr = multiply_lower_tri_self_transpose(L_corr);
    cov_matrix[J]  Sigma = quad_form_diag(corr, sigma);
    
    vector[J] obs_rep[N];
    vector[N] log_lik;
    
    if (get_obs_rep == 1) {
        for (n in 1:N) {
            obs_rep[n] = multi_normal_cholesky_rng(mu_cond[cond_id[n],sub_id[n]],
              diag_pre_multiply(sigma, L_corr));
        }    
    }
    if (get_log_lik == 1) {
        for (n in 1:N) {
            log_lik[n] = multi_normal_cholesky_lpdf(obs[n] | mu_cond[cond_id[n],sub_id[n]], 
              diag_pre_multiply(sigma, L_corr));
        }
    }
}
