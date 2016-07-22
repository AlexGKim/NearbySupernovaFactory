data {
  int D;                // Number of supernovae
  int N_mags;
  int N_EWs;
  vector[N_mags] mag_obs[D];
  vector[N_EWs] EW_obs[D];
  matrix[N_mags, N_mags] mag_cov[D];
  matrix[N_EWs, N_EWs] EW_cov[D];
}

parameters {

  vector<lower=-150, upper=150>[2] EW[D];
  vector<lower=-2, upper=2>[N_mags] mag_int[D];

  # vector<lower=-0.2, upper=0.2>[5] c;
  vector<lower=-0.01, upper=0.01>[5] alpha;
  vector<lower=-0.1, upper=0.1>[5] beta;

  vector<lower=0., upper=3>[4] gamma_;

  cholesky_factor_corr[N_mags] L_Omega;
  vector<lower=0.0, upper = 0.1>[N_mags] L_sigma;

  simplex[D] Delta_unit;
  real <lower = 0, upper = 3*D/4.> Delta_scale;

  simplex[D] k_unit;
  real <lower = 0, upper = 3*D/4.> k_scale;
}

transformed parameters {
  vector[D] Delta;
  vector[5] gamma;
  vector[D] k;

  gamma[1] <- gamma_[1];
  gamma[2] <- gamma_[2];
  gamma[3] <- 1.;
  gamma[4] <- gamma_[3];
  gamma[5] <- gamma_[4];

  Delta <- Delta_scale*(Delta_unit - 1./D);
  k <- k_scale*(k_unit - 1./D);
}

model {
  vector[5] means[D];
  matrix[5,5] L_Sigma;

  for (d in 1:D) {
      # means[d] <- Delta[d] + c + alpha*EW[d,1]  + beta*EW[d,2];
      means[d] <- Delta[d] + alpha*EW[d,1]  + beta*EW[d,2];
  }

  increment_log_prob(cauchy_log(L_sigma, 0.1,2.5));
  increment_log_prob(lkj_corr_cholesky_log(L_Omega, 2.));
  L_Sigma <- diag_pre_multiply(L_sigma, L_Omega);

  for (d in 1:D) {
    increment_log_prob(multi_normal_cholesky_log(mag_int[d], means[d], L_Sigma));
    # increment_log_prob(normal_log(mag_int[d], means[d], L_sigma));
    increment_log_prob(multi_normal_log(mag_obs[d],mag_int[d]+gamma*k[d], mag_cov[d]));
    increment_log_prob(multi_normal_log(EW_obs[d],EW[d], EW_cov[d]));
  }  
}