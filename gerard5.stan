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

  vector<lower=-1, upper=1>[5] c;
  vector<lower=-0.02, upper=0.02>[5] alpha;
  vector<lower=-0.05, upper=0.15>[5] beta;
  cholesky_factor_corr[N_mags] L_Omega;
  vector<lower=0.0, upper = 0.2>[N_mags] L_sigma;

  vector<lower=0.05, upper=5>[4] gamma;
  simplex [D] k_simplex;
  real <lower = 0, upper=D> k_scale;

  simplex [D] Delta;
  real <lower = 0, upper = D> Delta_scale;

  vector<lower=-10, upper=10>[4] rho0;
  vector<lower=-10, upper=10>[4] rho1;
  simplex[D] R;
  real <lower = 0, upper = 10*D> R_scale;

  real <lower=-5, upper=5> R_mn;
  # vector<lower=-5,upper=5> [D] Delta;
  # real <lower=0, upper=0.2> sigma_Delta;
}


model {
  vector[5] means[D];
  vector[5] gamma_;
  vector[5] rho0_;
  vector[5] rho1_;
  vector [D] k_;
  vector[D] Delta_;
  vector[D] R_;
  matrix[5,5] L_Sigma;

  gamma_[1]<- gamma[1];
  gamma_[2] <- gamma[2];
  gamma_[3] <- 1.;
  gamma_[4] <- gamma[3];
  gamma_[5] <- gamma[4];

  rho0_[1]<- rho0[1];
  rho0_[2] <- rho0[2];
  rho0_[3] <- rho0[3];
  rho0_[4] <- rho0[4];
  rho0_[5] <- 0;


  rho1_[1]<- rho1[1];
  rho1_[2] <- rho1[2];
  rho1_[3] <- rho1[3];
  rho1_[4] <- rho1[4];
  rho1_[5] <- 0;

  k_ <- k_scale*(k_simplex-1./D);
  Delta_ <- Delta_scale*(Delta - 1./D);
  R_ <- R_scale*(R - 1./D);

  # k_ <- k_ - mean(k_);

  for (d in 1:D) {
      means[d] <- Delta_[d] + c + alpha*EW[d,1]  + beta*EW[d,2] + rho0_*R_[d] + rho1_*pow(R_[d],2);
  }

  L_sigma ~ cauchy(0.1,2.5);
  L_Omega ~ lkj_corr_cholesky(20.);
  L_Sigma <- diag_pre_multiply(L_sigma, L_Omega);

  # Delta ~ normal(0,sigma_Delta);

  for (d in 1:D) {
    mag_int[d] ~ multi_normal_cholesky(means[d], L_Sigma);
    mag_obs[d] ~ multi_normal(mag_int[d]+gamma_*k_[d], mag_cov[d]);
    EW_obs[d] ~ multi_normal(EW[d], EW_cov[d]);
  }  

  R_~normal(R_mn,1);

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