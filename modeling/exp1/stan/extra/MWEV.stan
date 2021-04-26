// Title: Univariate Normal Non-linear Exepected Utility Maximization Model
// Author: Tyler Adkins
// Date: June 25th, 2020

functions {
  int get_max_idx(real[] EU);
  #include "/Users/adkinsty/Box/LeeLab/Experiments/Exp_files/reach/stan/hpp/get_max_idx.hpp"
  vector wp(vector p,real gamma) {
    // computes subjective weighted probabilities, given gamma
    // formula from Prelec article
    vector[rows(p)] wp;
    for (i in 1:rows(p)) {
      wp[i] = exp(pow(-log(p[i]),gamma));  
    }
    return wp;
  }
  // vector wp(vector p,real gamma) {
  //   // computes subjective weighted probabilities, given gamma
  //   // formula from Gonzales article
  //   vector[rows(p)] wp;
  //   for (i in 1:rows(p)) {
  //     wp[i] = (pow(p[i],gamma)/(pow(p[i],gamma)+pow((1-p[i]),gamma))^(1/gamma));
  //   }
  //   return wp;
  // }
  real get_mu(vector aim, vector pg, vector pl, vector pgl, real l, real gamma) {
    // Returns point that maximizes subjective weighted utility from a grid of aimpoints
    // subjective weigted utilities
    real SWU[rows(aim)];
    // subjective weighted probabiltiies
    vector[rows(aim)] wpg; 
    vector[rows(aim)] wpl; 
    vector[rows(aim)] wpgl; 
    
    // control to minimize computation
    if (l == 0) {
      wpg = wp(pg,gamma);
      wpgl = wp(pgl,gamma);
      SWU = to_array_1d(wpg*100 + wpgl*100); // wpg*1 + wpgl*1 + wpl*0 + wpn*0
    } else if (l == -1) {
      wpg = wp(pg,gamma);
      wpl = wp(pl,gamma);
      SWU = to_array_1d(wpg*100 - wpl*100); // wpg*1 + wpgl*0 + wpl*l + wpn*0
    } else {
      wpg = wp(pg,gamma);
      wpl = wp(pl,gamma);
      wpgl = wp(pgl,gamma);
      SWU = to_array_1d(wpg*100 - wpgl*400 - wpl*500); // wpg*1 + wpgl*(g+l) + wpl*l + wpn*0
    }
    return aim[get_max_idx(SWU)];
  }
}
data {
  int<lower=0> N; // number of observations
  int<lower=1> G; // number of aimpoints to consider (in simulation)
  int<lower=1> C;  // number of unique conditions (ie loss values)
  int<lower=1> J;  // number of unique subjects
  int<lower=1> L; // num of subject covariates (eg traits)
  int<lower=1> jj[N];  // subject IDs
  int<lower=1> cc[N];  // condition IDs (1,2,3)
  vector[C] cond; // list of unique conditions (ie loss values)
  vector[G] aim; // possible aim-points
  matrix<lower=0,upper=1>[G,J] pg; // prob of gain only
  matrix<lower=0,upper=1>[G,J] pl; // prob of loss loss only
  matrix<lower=0,upper=1>[G,J] pgl; // prob of gain and loss (overlap)
  vector[N] y; // observations  
  // generated quantities, y/n
  int<lower=0,upper=1> get_yrep;
  int<lower=0,upper=1> get_log_lik;
  int<lower=0,upper=1> prior_only;  
}
parameters {
  real<lower=0,upper=1> gamma_m;              // group mean sigma
  real<lower=0,upper=1> gamma_s;              // sd of sigma between subjects
  vector<lower=0,upper=1>[J] gamma;
  real<lower=0> sigma_m;              // group mean sigma
  real<lower=0> sigma_s;              // sd of sigma between subjects
  vector<lower=0>[J] sigma;
}
transformed parameters {
  matrix[C, J] mu_cond;  // condition/subject specific aimpoints
  for (j in 1:J)
    for (c in 1:C)
      mu_cond[c,j] = get_mu(aim, pg[,j], pl[,j], pgl[,j], cond[c], gamma[j]);
}
model {
  gamma_m ~ beta(1,1);
  gamma_s ~ beta(1,1);
  gamma ~ normal(gamma_m, gamma_s);
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
