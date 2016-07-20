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

  vector<lower=-500, upper=500>[2] EW[D];
  vector<lower=-5, upper=5>[N_mags] mag_int[D];

  vector<lower=-0.2, upper=0.2>[5] c;
  vector<lower=-0.02, upper=0.02>[5] alpha;
  vector<lower=-0.2, upper=0.2>[5] beta;

  vector<lower=0, upper=10>[4] gamma_;

  cholesky_factor_corr[N_mags] L_Omega;
  vector<lower=0.0, upper = 0.15>[N_mags] L_sigma;

  simplex[D] Delta_unit;
  # simplex [D] Delta_simplex;
  real <lower = 0, upper = D/2.> Delta_scale;

  simplex[D] k_unit;
  real <lower = 0, upper = D/2.> k_scale;
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
  # for (d in 1:D-1) {
  #   k[d] <- k_unit[d];
  #   Delta[d] <- Delta_unit[d];
  # }
  # k[D] <- -sum(k_unit);
  # Delta[D] <- -sum(Delta_unit);

  # k <- rotation * k;
  # Delta <- rotation * Delta;
}

model {
  vector[5] means[D];
  matrix[5,5] L_Sigma;

  for (d in 1:D) {
      means[d] <- Delta[d] + c + alpha*EW[d,1]  + beta*EW[d,2];
  }

  L_sigma ~ cauchy(0.1,2.5);
  L_Omega ~ lkj_corr_cholesky(2.);
  L_Sigma <- diag_pre_multiply(L_sigma, L_Omega);


  for (d in 1:D) {
    mag_int[d] ~ multi_normal_cholesky(means[d], L_Sigma);
    # mag_int[d] ~ normal(means[d], L_sigma);

    mag_obs[d] ~ multi_normal(mag_int[d]+gamma*k[d], mag_cov[d]);
    EW_obs[d] ~ multi_normal(EW[d], EW_cov[d]);
  }  
}