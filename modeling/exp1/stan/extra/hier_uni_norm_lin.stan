// Title: Hierarchical Univariate Normal Linear Regression
// Author: Tyler Adkins
// Date: June 19th, 2020

data {
    int<lower=0> N; // number of observations
    int<lower=2,upper=2> J;  // number of variates
    vector[N] ratio; // gain covariate
    row_vector[J] obs[N]; // observations
    int<lower=0,upper=1> get_obs_rep;  // get posterior predictions?
    int<lower=0,upper=1> get_log_lik;  // get log liklihoods?
    int<lower=0,upper=1> prior_only;   // ignore the data?
}
parameters {
    row_vector[J] alpha;  // intercepts
    row_vector[J] beta;   // slopes
    vector<lower=0>[J] sigma;
    cholesky_factor_corr[J] L_corr;
}

model {
    row_vector[J] mu[N];
    // prior
    alpha ~ normal(0,1); 
    beta ~ normal(0,1);
    sigma ~ normal(0, 1);
    L_corr ~ lkj_corr_cholesky(2.0);
    // likelihood
    if (prior_only == 0) {
        for (n in 1:N) {
            mu[n] = alpha + ratio[n]*beta;
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
                alpha + ratio[n]*beta, diag_pre_multiply(sigma, L_corr));
        }    
    }
    if (get_log_lik == 1) {
        for (n in 1:N) {
            log_lik[n] = multi_normal_cholesky_lpdf(
                obs[n] | alpha + ratio[n]*beta, diag_pre_multiply(sigma, L_corr));
        }
    }

}

