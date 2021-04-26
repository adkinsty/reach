// Title: Univariate Normal Non-linear Exepected Utility Maximization Model
// Author: Tyler Adkins
// Date: June 25th, 2020

functions {
  int get_max_idx(real[] EU);
  #include "/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/stan/hpp/get_max_idx.hpp"
  real get_mu(vector aim, vector pg, vector pl, vector pgl, real l, real alpha, real beta) {
    // Get aimpoint that maximizes expected utility
    real EU[rows(aim)];
    real ug;
    real ul;
    real ugl;
    // scale payoffs by 100x to avoid multiplication by 1
    if (l==0) {
      ug = pow(100,alpha);
      EU = to_array_1d(pg*ug + pgl*ug); // + pl*0 + pn*0
    } else if (l==-1) {
      ug = pow(100,alpha);
      ul = -pow(100,beta);
      EU = to_array_1d(pg*ug + pl*ul); // + pgl*0 + pn*0
    } else {
      ug = pow(100,alpha);
      ul = -pow(500,beta);
      ugl = -pow(400,beta);
      EU = to_array_1d(pg*ug + pgl*ugl + pl*ul); // + pn*0
    }
    return aim[get_max_idx(EU)];
  }
}
data {
  int<lower=0> N; // number of observations
  int<lower=1> G; // number of aimpoints to consider
  int<lower=1> C;  // number of unique conditions
  int<lower=1> J;  // number of unique subjects
  int<lower=1> jj[N];  // subject IDs
  int<lower=1> cc[N];  // condition IDs (1,2,3)
  vector[C] cond; // list of unique conditions
  vector[G] aim; // possible aim-points
  matrix<lower=0,upper=1>[G,J] pg; // prob of gain
  matrix<lower=0,upper=1>[G,J] pl; // prob of loss
  matrix<lower=0,upper=1>[G,J] pgl; // prob of gain and loss
  vector[N] y; // observations  
  // generated quantities, y/n
  int<lower=0,upper=1> get_yrep;
  int<lower=0,upper=1> get_log_lik;
  int<lower=0,upper=1> prior_only;  
}
parameters {
  real<lower=0> sigma_m;              // group mean sigma
  real<lower=0> sigma_s;              // sd of sigma between subjects
  vector<lower=0>[J] sigma;
}
transformed parameters {
  matrix[C, J] mu_cond;  // condition/subject specific aimpoints
  for (j in 1:J)
    for (c in 1:C)
      mu_cond[c,j] = get_mu(aim, pg[,j], pl[,j], pgl[,j], cond[c]);
}
model {
  sigma_m ~ std_normal();
  sigma_s ~ std_normal();
  sigma ~ normal(sigma_m,sigma_s);

  if (prior_only == 0) {
    {
      vector[N] Mu;
      vector[N] Sigma;
      for (n in 1:N) {
        Mu[n] = mu_cond[cc[n],jj[n]];
        Sigma[n] = sigma[jj[n]];
      }
      y ~ normal(Mu,Sigma); // likelihood
    }
  }
}
generated quantities {
    vector[N] yrep;
    vector[N] log_lik;
    
    if (get_yrep == 1)
        for (n in 1:N)
            yrep[n] = normal_rng(mu_cond[cc[n],jj[n]],sigma[jj[n]]);
    if (get_log_lik == 1)
        for (n in 1:N)
            log_lik[n] = normal_lpdf(y[n] | mu_cond[cc[n],jj[n]],sigma[jj[n]]);
}
