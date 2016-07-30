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

  vector<lower=-.1, upper=.1>[5] c;
  vector<lower=-0.01, upper=0.01>[5] alpha;
  vector<lower=-0.1, upper=0.1>[5] beta;

  real<lower=0,upper=.5> gamma0;
  vector<lower=-0.5, upper=0.5>[4] gamma_;

  simplex[D] k_unit;

  cholesky_factor_corr[N_mags] L_Omega;
  vector<lower=0.0, upper = 0.12>[N_mags] L_sigma;

  simplex[D] Delta_unit;
  real <lower = 0, upper = D/4.> Delta_scale;

  simplex[D] R_unit;
  # real<lower=0, upper=0.5> rho00;
  # vector<lower=-0.4, upper=0.5>[4] rho0_;
  vector<lower=-0.05, upper=0.05>[5] rho1;
}

transformed parameters {
  vector[D] Delta;
  vector[5] gamma;
  # vector[5] rho0;
  vector[D] R;
  vector[D] k;

  gamma[1] <- gamma0;
  gamma[2] <- gamma_[1];
  gamma[3] <- gamma_[2];
  gamma[4] <- gamma_[3];
  gamma[5] <- gamma_[4];

  # rho0[1] <- rho0_[1];
  # rho0[2] <- rho0_[2];
  # rho0[3] <- rho0_[3];
  # rho0[4] <- rho0_[4];
  # rho0[5] <- rho00;


  Delta <- Delta_scale*(Delta_unit - 1./D);
  k <- (k_unit - 1./D);
  k <- k / sqrt(sum(k .* k)/D);
  R <- (R_unit - mean(R_unit));
  R <- R_unit / sqrt(sum(R .* R)/D);
}

model {
  vector[5] means[D];
  matrix[5,5] L_Sigma;

  for (d in 1:D) {
      means[d] <- Delta[d] + c+ alpha*EW[d,1] + beta*EW[d,2] + rho1*(pow(R[d],2));
  }

  increment_log_prob(cauchy_log(L_sigma, 0.,1));
  increment_log_prob(lkj_corr_cholesky_log(L_Omega, 4.));
  L_Sigma <- diag_pre_multiply(L_sigma, L_Omega);

  for (d in 1:D) {

    increment_log_prob(multi_normal_cholesky_log(mag_int[d], means[d], L_Sigma));
    # increment_log_prob(normal_log(mag_int[d], means[d], L_sigma));
    increment_log_prob(multi_normal_log(mag_obs[d],mag_int[d]+gamma*k[d], mag_cov[d]));
    increment_log_prob(multi_normal_log(EW_obs[d],EW[d], EW_cov[d]));

  }  
}