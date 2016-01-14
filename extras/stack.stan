// Code based on Version 2.9.0 of the Stan manual
// Section 6.14. Multivariate Outcomes, Multivariate Probit Regression
// Licensed under CC-BY (https://creativecommons.org/licenses/by/4.0/legalcode)
// Note that most of the data manipulation code is not included here because
// operations like summing over integer arrays are easier in R than Stan.
data {
  // Dimensions of training data
  int<lower=1> D;                       // number of columns (species)
  int<lower=0> N;                       // number of rows (sites)

  // Presence/absence array
  int<lower=0,upper=1> y[N,D];          // binary presence/absence array

  // Presence/absence sums
  int<lower=0,upper=N*D> N_pos;         // numer of presences in y
  int<lower=0,upper=N*D> N_neg;         // numer of absences in y

  // row and column indices for y==1
  int<lower=1,upper=N> n_pos[N_pos];    // row(y)[y==1]
  int<lower=1,upper=D> d_pos[N_pos];    // col(y)[y==1]

  // row and column indices for y==0
  int<lower=1,upper=N> n_neg[N_neg];    // row(y)[y==0]
  int<lower=1,upper=D> d_neg[N_neg];    // col(y)[y==0]

  // Import the SDM predictions
  vector[D] probit_p[N];  // probit(Predicted occurrence probabilities)

  // Priors
  real<lower=0> lkj_eta; // lkj prior shape parameter for correlation strength
}

transformed data {
  vector[D] scaled_probit_p[N];

  // Fix discrepancy between p and expected value of
  // inv_probit(probit_p + rnorm(length(p))).
  // Dividing by sqrt(2) works empirically, but will need a clearer explanation
  // for what's going on in this step.
  for (i in 1:N) {
    scaled_probit_p[i] <- probit_p[i] / sqrt(2);
  }
}

parameters{
  cholesky_factor_corr[D] L_rho;
  vector<lower=0>[N_pos] z_pos; // Sampled values for latent Gaussian given y==1
  vector<upper=0>[N_neg] z_neg; // Sampled values for latent Gaussian given y==0
}

transformed parameters {
  vector[D] z[N];               // mean of the latent normal

  // Combine z_pos and z_neg into a single array called z
  for (n in 1:N_pos){
    z[n_pos[n], d_pos[n]] <- z_pos[n];
  }
  for (n in 1:N_neg){
    z[n_neg[n], d_neg[n]] <- z_neg[n];
  }
}

model {
  // Sample a factorized correlation matrix from posterior
  L_rho ~ lkj_corr_cholesky(lkj_eta);

  // Sample Gaussian random deviates from posterior
  for (i in 1:N){
    z[i] ~ multi_normal_cholesky(scaled_probit_p[i], L_rho);
  }
}

generated quantities {
  corr_matrix[D] rho;

  // Undo Cholesky factorization
  rho <- multiply_lower_tri_self_transpose(L_rho);
}
