data {
  int D;                // Number of supernovae
  int N_colors;
  int N_EWs;
  vector[N_colors] color_obs[D];
  vector[N_EWs] EW_obs[D];
  matrix[N_colors, N_colors] color_cov[D];
  matrix[N_EWs, N_EWs] EW_cov[D];

  real color_std;
  vector[N_colors] color_mn;
}

parameters {

  vector<lower=-5, upper=5>[2] EW[D];

  vector[N_colors] color_int[D];
  vector<lower=-1, upper=1>[5] c;
  vector<lower=-5, upper=5>[5] alpha;
  vector<lower=-5, upper=5>[5] beta;
  cholesky_factor_corr[N_colors] L_Omega;
  vector<lower=0>[N_colors] L_sigma;

  vector<lower=-5, upper=5>[4] gamma;
  real<lower=-5, upper=5> k[D-1];
}

model {
  vector[5] means[D];
  vector[5] gamma_;
  real k_[D];
  matrix[5,5] L_Sigma;

  gamma_[1]<- gamma[1];
  gamma_[2] <- 1.;
  gamma_[3] <- gamma[2];
  gamma_[4] <- gamma[3];
  gamma_[5] <- gamma[4];

  k_[1] <- 0;
  for (d in 1:(D-2)) {
    k_[1+d] <- k[d];
  }

  for (d in 1:(D-1)) {
      means[d] <- c + alpha*EW[d,1]  + beta*EW[d,2];
  }

  L_sigma ~ cauchy(0,2.5);
  L_Omega ~ lkj_corr_cholesky(2.);
  L_Sigma <- diag_pre_multiply(L_sigma, L_Omega);

  for (d in 1:(D-1)) {
    color_int[d] ~ multi_normal_cholesky(means[d], L_Sigma);
    color_obs[d] ~ multi_normal(means[d], color_cov[d]);
    EW_obs[d] ~ multi_normal(EW[d], EW_cov[d]);
  }  
}

generated quantities {
  # corr_matrix[4] Omega_color;
  vector[5] c_unnorm;
  vector[5] alpha_unnorm;
  vector[5] beta_unnorm;
  vector[4] gamma_unnorm;
  # vector[3] gamma_unnorm;
  c_unnorm <- c*color_std + color_mn;
  alpha_unnorm <- alpha*color_std;
  beta_unnorm <- beta*color_std;
  gamma_unnorm <- gamma*color_std;
}
