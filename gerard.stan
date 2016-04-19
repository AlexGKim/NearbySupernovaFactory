data {
  int D;                // Number of supernovae
  int N_colors;
  int N_EWs;
  vector[N_colors] color_obs[D];
  vector[N_EWs] EW_obs[D];
  matrix[N_colors, N_colors] color_cov[D];
  matrix[N_EWs, N_EWs] EW_cov[D];
}

parameters {
  vector<lower=-2, upper=2>[2] EW0;
  vector<lower=0, upper=2>[2] L_sigma_EW;
  cholesky_factor_corr[2] L_Omega_EW;
  vector<lower=-5, upper=5>[2] EW[D];

  vector<lower=-2, upper=2>[4] c;
  vector<lower=-10, upper=10>[4] alpha;
  vector<lower=-10, upper=10>[4] beta;
  vector<lower=-10, upper=10>[3] gamma;

  vector<lower=0, upper=2>[4] L_sigma_color;
  cholesky_factor_corr[4] L_Omega_color;
  real<lower=-5, upper=5> k[D];
  vector<lower=0, upper=5>[4] EXY[D];
  vector<lower=-5, upper=5>[4] colors[D];
}

transformed parameters{
  vector[4] means[D];
  vector[4] gamma_;
  gamma_[1]<- gamma[1];
  gamma_[2] <- 1.;
  gamma_[3] <- gamma[2];
  gamma_[4] <- gamma[3];
  for (d in 1:(D-1)) {
      means[d] <- c + alpha*EW[d,1]  + beta*EW[d,2] + gamma_*k[d] + colors[d];
  }
}

model {
  matrix[2,2] L_Sigma_EW;
  matrix[4,4] L_Sigma_color;

  L_sigma_EW ~ cauchy(0,2.5);
  L_Omega_EW ~ lkj_corr_cholesky(2.);
  L_sigma_color ~ cauchy(0,2.5);
  L_Omega_color ~ lkj_corr_cholesky(2.);
  L_Sigma_EW <- diag_pre_multiply(L_sigma_EW, L_Omega_EW);
  L_Sigma_color <- diag_pre_multiply(L_sigma_color, L_Omega_color);

  for (d in 1:(D-1)) {
    EW[d] ~ multi_normal_cholesky(EW0,L_Sigma_EW);
    colors[d] ~ multi_normal_cholesky(means[d], L_Sigma_color);
    color_obs[d] ~ multi_normal(colors[d], color_cov[d]);
    EW_obs[d] ~ multi_normal(EW[d], EW_cov[d]);
  }  
}