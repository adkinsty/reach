// Title: Perceptual heuristic model
// Author: Tyler Adkins
// Date: June 25th, 2020

data {
  int<lower=0> N; // number of observations
  int<lower=1> C;  // number of unique conditions
  int<lower=1> J;  // number of unique subjects
  int<lower=1> jj[N];  // subject IDs
  int<lower=1> cc[N];  // condition IDs (1,2,3)
  matrix[C,J] aim;   // optimal aimpoints
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
model {
  // priors
  sigma_m ~ std_normal();
  sigma_s ~ std_normal();
  sigma ~ normal(sigma_m, sigma_s);
  
  if (prior_only == 0) {
    {
      vector[N] Mu;
      vector[N] Sigma;
      // set up predictors
      for (n in 1:N) {
        Mu[n] = aim[cc[n],jj[n]];
        Sigma[n] = sigma[jj[n]];
      }
      // get likelihood
      y ~ normal(Mu,Sigma);
    }
  }
}
generated quantities {
    vector[N] yrep;
    vector[N] log_lik;
    
    if (get_yrep == 1) {
        for (n in 1:N) {
            yrep[n] = normal_rng(aim[cc[n],jj[n]],sigma[jj[n]]);
        }    
    }
    if (get_log_lik == 1) {
        for (n in 1:N) {
            log_lik[n] = normal_lpdf(y[n] | aim[cc[n],jj[n]],sigma[jj[n]]);
        }
    }
}
