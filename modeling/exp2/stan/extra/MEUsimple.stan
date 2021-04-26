// Title: Univariate Normal Non-linear Exepected Utility Maximization Model
// Author: Tyler Adkins
// Date: June 25th, 2020

functions {
  int get_max_idx(real[] EU);
  #include "/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/stan/hpp/get_max_idx.hpp"
  vector wp(vector p,real gamma) {
    vector[rows(p)] wp;
    for (i in 1:rows(p))
      wp[i] = (pow(p[i],gamma)/(pow(p[i],gamma)+pow((1-p[i]),gamma))^(1/gamma));
    return wp;
  }
  real get_mu(vector aim, vector pg, vector pl, vector pgl, vector pn, real g, real l, real beta, real gamma) {
    // Get aimpoint that maximizes expected utility
    real EU[rows(aim)];
    real ul;
    real wpg; 
    real wpl; 
    real wpgl; 
    real wpn;
    wpg = wp(pg,gamma);
    wpl = wp(pl,gamma);
    wpgl = wp(pgl,gamma);
    wpn = wp(pgl,gamma);
    ul = -pow(-l,beta);
    EU = to_array_1d(wpg*g + wpgl*(g+ul) + wpl*ul + wpn*0);
    return aim[get_max_idx(EU)];
  }
}
data {
  int<lower=0> N; // number of observations
  int<lower=1> G; // number of aimpoints to consider
  int<lower=1> C;  // number of unique conditions
  int<lower=1> J;  // number of unique subjects
  int<lower=1> L; // num of subject covariates
  int<lower=1> jj[N];  // subject IDs
  int<lower=1> cc[N];  // condition IDs (1,2,3)
  vector[C] cond; // list of unique conditions
  vector[N] gain; // gain covariate
  vector[N] loss; // gain covariate
  vector[G] aim; // possible aim-points
  matrix<lower=0,upper=1>[G,J] pg; // prob of gain
  matrix<lower=0,upper=1>[G,J] pl; // prob of loss
  matrix<lower=0,upper=1>[G,J] pgl; // prob of gain and loss
  matrix<lower=0,upper=1>[G,J] pn; // prob of neither
  vector[N] y; // observations  
  // generated quantities, y/n
  int<lower=0,upper=1> get_yrep;
  int<lower=0,upper=1> get_log_lik;
  int<lower=0,upper=1> prior_only;  
}
parameters {
  real<lower=0> beta_m;              // group mean sigma
  real<lower=0> beta_s;              // sd of sigma between subjects
  vector<lower=0.000001>[J] beta;
  real<lower=0> gamma_m;              // group mean sigma
  real<lower=0> gamma_s;              // sd of sigma between subjects
  vector<lower=0.000001>[J] gamma;
  real<lower=0> sigma_m;              // group mean sigma
  real<lower=0> sigma_s;              // sd of sigma between subjects
  vector<lower=0>[J] sigma;
}
transformed parameters {
  matrix[C, J] mu_cond;  // condition/subject specific aimpoints
  for (j in 1:J)
    for (c in 1:C)
      mu_cond[c,j] = get_mu(aim, pg[,j], pl[,j], pgl[,j], pn[,j], 1, cond[c], beta[j]);
}
model {
  beta_m ~ gamma(10,10);
  beta_s ~ std_normal();
  beta ~ normal(beta_m, beta_s);
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
