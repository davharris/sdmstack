// Code based on Version 2.9.0 of the Stan manual
// Section 6.14. Multivariate Outcomes, Multivariate Probit Regression
// Distributed under CC-BY  license (https://creativecommons.org/licenses/by/4.0/legalcode)
// Note that most of the data manipulation code is not included here because 
// operations like summing over integer arrays are easier in R than Stan.
data {
  // Dimensions
  int<lower=1> D;                       // number of columns (species)
  int<lower=0> N;                       // number of rows (sites)
  int<lower=0,upper=1> y[N,D];          // binary presence/absence array
  int<lower=0,upper=N*D> N_pos;         // numer of presences in y
  int<lower=0,upper=N*D> N_neg;         // numer of absences in y
  
  // row and column indices for y==1
  int<lower=1,upper=N> n_pos[N_pos];    // rows
  int<lower=1,upper=D> d_pos[N_pos];    // columns
  
  // row and column indices for y==0
  int<lower=1,upper=N> n_neg[N_neg];    // rows
  int<lower=1,upper=D> d_neg[N_neg];    // columns
  
  // Import the SDM predictions
  vector[D] probit_p[N];  // probit(Predicted occurrence probabilities)
  
  // Priors
  real<lower=0> lkj_eta; // lkj prior shape parameter for correlation strength
}

parameters{
  cholesky_factor_corr[D] L_rho;
  vector<lower=0>[N_pos] z_pos; // Sampled values for latent Gaussian for y==1
  vector<upper=0>[N_neg] z_neg; // Sampled values for latent Gaussian for y==0
  
  vector[D] intercepts;
}

transformed parameters {
  vector[D] z[N]; 
  vector[D] z_minus_predictions[N]; 
  
  // Combine z_pos and z_neg into a single vector called z
  for (n in 1:N_pos){
    z[n_pos[n], d_pos[n]] <- z_pos[n];
  }
  for (n in 1:N_neg){
    z[n_neg[n], d_neg[n]] <- z_neg[n];
  }
  
  // full z is z_minus_predictions plus the predictions from the SDM models
  for (i in 1:N){
    z_minus_predictions[i] <- z[i] - probit_p[i];
  }
}

model {
  intercepts ~ normal(0.0, 5);
  L_rho ~ lkj_corr_cholesky(lkj_eta);
  z_minus_predictions ~ multi_normal_cholesky(intercepts, L_rho);
}

generated quantities {
  corr_matrix[D] rho;
  rho <- multiply_lower_tri_self_transpose(L_rho);
}
