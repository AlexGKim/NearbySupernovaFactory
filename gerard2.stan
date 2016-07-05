data {
  int D;                // Number of supernovae
  int N_mags;
  int N_EWs;
  vector[N_mags] mag_obs[D];
  vector[N_EWs] EW_obs[D];
  matrix[N_mags, N_mags] mag_cov[D];
  matrix[N_EWs, N_EWs] EW_cov[D];

  # real mag_std;
  # vector[N_mags] mag_mn;
}

parameters {

  vector<lower=-100, upper=100>[2] EW[D];
  vector<lower=-2, upper=2>[N_mags] mag_int[D];

  vector<lower=-0.3, upper=0.3>[5] c;
  vector<lower=-0.05, upper=0.05>[5] alpha;
  vector<lower=0.005, upper=0.1>[5] beta;
  cholesky_factor_corr[N_mags] L_Omega;
  vector<lower=0.02, upper = 0.2>[N_mags] L_sigma;

  vector<lower=0.5, upper=1.8>[4] gamma;
  real<lower=-3, upper=3> k[D-1];
}

model {
  vector[5] means[D];
  vector[5] gamma_;
  real k_[D];
  matrix[5,5] L_Sigma;

  gamma_[1]<- gamma[1];
  gamma_[2] <- gamma[2];
  gamma_[3] <- 1.;
  gamma_[4] <- gamma[3];
  gamma_[5] <- gamma[4];

  k_[1] <- 0;
  for (d in 1:(D-2)) {
    k_[1+d] <- k[d];
  }

  for (d in 1:(D-1)) {
      means[d] <- c + alpha*EW[d,1]  + beta*EW[d,2];
  }

  L_sigma ~ cauchy(0.1,2.5);
  L_Omega ~ lkj_corr_cholesky(2.);
  L_Sigma <- diag_pre_multiply(L_sigma, L_Omega);

  for (d in 1:(D-1)) {
    mag_int[d] ~ multi_normal_cholesky(means[d], L_Sigma);
    mag_obs[d] ~ multi_normal(mag_int[d]+gamma_*k[d], mag_cov[d]);
    EW_obs[d] ~ multi_normal(EW[d], EW_cov[d]);
  }  
}

# generated quantities {
#   # corr_matrix[4] Omega_mags;
#   vector[5] c_unnorm;
#   vector[5] alpha_unnorm;
#   vector[5] beta_unnorm;
#   vector[4] gamma_unnorm;
#   # vector[3] gamma_unnorm;
#   c_unnorm <- c*mags_std + mags_mn;
#   alpha_unnorm <- alpha*mags_std;
#   beta_unnorm <- beta*mags_std;
#   gamma_unnorm <- gamma*mags_std;
# }
