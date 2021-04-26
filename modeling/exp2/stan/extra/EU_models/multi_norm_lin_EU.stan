// Title: Maximum Exepected Utility Model
// Author: Tyler Adkins
// Date: June 19th, 2020


functions {
  real u(real v, real alpha, real beta, real lambda){
    // Utility function
    if (v >= 0) {
      return v^alpha;
    } else {
      return -lambda*(-v)^beta;
    }
  }
  vector get_aim(vector grid, vector p_gain, vector p_loss, real gain, 
                 real loss, real alpha, real beta, real lambda) {
    // Get aimpoint that maximizes expected utility
    vector[2] aim;
    vector[rows(grid)] EU;

    EU = p_gain * u(gain,alpha,beta,lambda) + p_loss * u(loss,alpha,beta,lambda);

    aim[1] = grid[sort_indices_asc(EU)[rows(grid)]]; // maximizing x aimpoint
    aim[2] = 0; // never aim above midline
    return aim;
  }
}

data {
  int<lower=0> N; // number of observations
  int<lower=1> G; // number of aimpoints to consider
  int<lower=2,upper=2> J;  // number of variates (x,y)
  int<lower=1> C; // number of unique conditions
  vector[C] cond; // list of unique conditions
  vector[N] gain; // gain covariate
  vector[N] loss; // gain covariate
  vector[G] grid; // candidate aim-points
  vector[G] p_gain; // gain probabilities for each aim-point
  vector[G] p_loss; // loss probabilities for each aim-point
  row_vector[J] obs[N]; // observations  // generated quantities, y/n
  int<lower=0,upper=1> get_obs_rep;
  int<lower=0,upper=1> get_log_lik;
  int<lower=0,upper=1> prior_only;  
}
parameters {
  real<lower=0> alpha;  // gain growth parameter
  real<lower=0> beta;   // loss growth parameter
  real<lower=0> lambda; // loss aversion parameter
  vector<lower=0>[J] sigma;
  cholesky_factor_corr[J] L_corr;
}
transformed parameters {
  vector[J] mu_cond[C];
  
  for (c in 1:C) {
    mu_cond[c] = get_aim(grid, p_gain, p_loss, 100, cond[c], alpha, beta, lambda);
  }
}
model {
  vector[J] mu[N];

  // priors
  alpha ~ normal(1,.5); 
  beta ~ normal(1,.5);
  lambda ~ normal(1,.5);
  sigma ~ normal(0, 1);
  L_corr ~ lkj_corr_cholesky(2.0);
  // likelihood
  if (prior_only == 0) {
    for (n in 1:N) {
      mu[n] = get_aim(grid, p_gain, p_loss, gain[n], loss[n], alpha, beta, lambda);
    }
    obs ~ multi_normal_cholesky(mu, diag_pre_multiply(sigma, L_corr));
  }
}
generated quantities {
    corr_matrix[J] corr = multiply_lower_tri_self_transpose(L_corr);
    cov_matrix[J]  Sigma = quad_form_diag(corr, sigma);
    
    vector[J] obs_rep[N];
    vector[N] log_lik;
    
    if (get_obs_rep == 1) {
        for (n in 1:N) {
            obs_rep[n] = multi_normal_cholesky_rng(
                get_aim(grid, p_gain, p_loss, gain[n], loss[n], alpha, beta, lambda),
                diag_pre_multiply(sigma, L_corr));
        }    
    }
    if (get_log_lik == 1) {
        for (n in 1:N) {
            log_lik[n] = multi_normal_cholesky_lpdf(
                obs[n] | get_aim(grid, p_gain, p_loss, gain[n], loss[n], alpha, beta, lambda),
                         diag_pre_multiply(sigma, L_corr));
        }
    }
}
