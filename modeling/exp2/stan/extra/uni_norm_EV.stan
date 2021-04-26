// Title: Univariate Maximum Exepected Value Model
// Author: Tyler Adkins
// Date: June 25th, 2020


functions {
  int get_max_idx(real[] EU);
  #include "/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/stan/hpp/get_max_idx.hpp"
  real get_aim(vector grid, vector p_gain, vector p_loss, int gain, real loss) {
    // Get aimpoint that maximizes expected utility
    real EU[rows(grid)];
    EU = to_array_1d(p_gain * gain + p_loss * loss);
    return grid[get_max_idx(EU)];
  }
}

data {
  int<lower=0> N; // number of observations
  int<lower=1> G; // number of aimpoints to consider
  int<lower=1> C;  // number of unique conditions
  int<lower=1> J;  // number of unique subjects
  int<lower=1> jj[N];  // subject IDs
  int<lower=1> cc[N];  // condition IDs (1,2,3)
  vector[C] cond; // list of unique losses, ordered from smallest to largest
  vector[N] gain; // gain covariate
  vector[N] loss; // gain covariate
  vector[G] grid; // candidate aim-points
  matrix<lower=0,upper=1>[G,J] p_gain; // gain probabilities for each aim-point
  matrix<lower=0,upper=1>[G,J] p_loss; // loss probabilities for each aim-point
  vector[N] y; // observations  
  // generated quantities, y/n
  int<lower=0,upper=1> get_yrep;
  int<lower=0,upper=1> get_log_lik;
  int<lower=0,upper=1> prior_only;  
}
parameters {
  real<lower=0> sigma_m;
  real<lower=0> sigma_s;
  vector<lower=0>[J] sigma;
}
transformed parameters {
  matrix[C,J] Mu;
  // pre-compute condition/subject means
  for (c in 1:C) {
    for (j in 1:J) {
      Mu[c,j] = get_aim(grid, p_gain[,j], p_loss[,j], 1, cond[c]);
    }
  }
}
model {
  // priors
  sigma_m ~ std_normal();
  sigma_s ~ std_normal();
  sigma ~ normal(sigma_m, sigma_s);
  
  if (prior_only == 0) {
    {
      vector[N] mu;
      vector[N] Sigma;
      // set up predictors
      for (n in 1:N) {
        mu[n] = Mu[cc[n],jj[n]];
      }
      // get likelihood
      y ~ normal(mu,sigma);
    }
  }
}
generated quantities {
    vector[N] yrep;
    vector[N] log_lik;
    
    if (get_yrep == 1) {
        for (n in 1:N) {
            yrep[n] = normal_rng(Mu[cc[n],jj[n]],sigma[jj[n]]);
        }    
    }
    if (get_log_lik == 1) {
        for (n in 1:N) {
            log_lik[n] = normal_lpdf(y[n] | Mu[cc[n],jj[n]],sigma[jj[n]]);
        }
    }
}
