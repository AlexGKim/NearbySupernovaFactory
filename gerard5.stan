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
  vector[N_mags] mag_int[D];

  vector[5] c;
  vector<lower=-0.02, upper=0.02>[5] alpha;
  vector<lower=-0.05, upper=0.15>[5] beta;

  real<lower=0> gamma0;
  vector[4] gamma_;
  real<lower=0> rho00;
  vector[4] rho0_;
  
  vector[5] rho1;

  cholesky_factor_corr[N_mags] L_Omega;
  vector<lower=0.0, upper = 0.5>[N_mags] L_sigma;

  simplex [D] Delta_simplex;
  real <lower = 0, upper = D> Delta_scale;

  vector<lower=-10, upper=10> [D-2] k_;
  vector<lower=-10, upper=10> [D-2] R_;
  # real <lower = 0, upper=D> k_scale;


  # real<lower=0, upper=5> rho00;


  # real<lower = 0, upper = D> R_scale;

  # real <lower=-5, upper=5> R_mn;
  # vector<lower=-5,upper=5> [D] Delta;
  # real <lower=0, upper=0.2> sigma_Delta;
}

transformed parameters {
  vector[D] Delta;
  vector[D] k;
  vector[D] R;

  vector[5] gamma;
  vector[5] rho0;

  rho0[1] <- rho00;
  rho0[2] <- rho0_[1];
  rho0[3] <- rho0_[2];
  rho0[4] <- rho0_[3];
  rho0[5] <- rho0_[4];

  gamma[1] <- gamma0;
  gamma[2] <- gamma_[1];
  gamma[3] <- gamma_[2];
  gamma[4] <- gamma_[3];
  gamma[5] <- gamma_[4];


  Delta <- Delta_scale*(Delta_simplex - 1./D);

  for (d in 1:D-2){
    k[d] <- k_[d];
    R[d] <- R_[d];
  }

  {
    real mean_;
    real variance_;
    real term;
    mean_ <- sum(k_);
    variance_ <- sum(k_ .* k_);
    term <- 0.5*sqrt(2*D + 2*variance_ - mean_*mean_);
    k[D-1] <- -0.5*mean_ + term;
    k[D] <- -0.5*mean_ - term;

    mean_ <- sum(R_);
    variance_ <- sum(R_ .* R_);
    term <- 0.5*sqrt(2*D + 2*variance_ - mean_*mean_);
    R[D-1] <- -0.5*mean_ + term;
    R[D] <- -0.5*mean_ - term;
  }
  # k <- (k_simplex-1./D);
  # k <- k/sqrt(sum(k .* k)/D);
  # R <- (R_simplex-1./D);
  # R <- R/sqrt(sum(R .* R)*D);
}

model {
  vector[5] means[D];
  # vector[5] gamma_;

  # vector[D] Delta_;
  # vector[D] R_;
  matrix[5,5] L_Sigma;

  # gamma_[1]<- gamma[1];
  # gamma_[2] <- gamma[2];
  # gamma_[3] <- 1.;
  # gamma_[4] <- gamma[3];
  # gamma_[5] <- gamma[4];

  # rho0_[1] <- rho0[1];
  # rho0_[2] <- rho0[2];
  # rho0_[3] <- 1;
  # rho0_[4] <- rho0[3];
  # rho0_[5] <- rho0[4];


  # rho1_[1]<- rho1[1];
  # rho1_[2] <- rho1[2];
  # rho1_[3] <- 1;
  # rho1_[4] <- rho1[3];
  # rho1_[5] <- rho1[4];


  # R_ <- R_scale*(R - 1./D);
  # R_ <- D * R - 1.;

  # k_ <- k_ - mean(k_);

  for (d in 1:D) {
      means[d] <- Delta[d] + c + alpha*EW[d,1]  + beta*EW[d,2] + rho0*R[d] + rho1*pow(R[d],2);
  }

  L_sigma ~ cauchy(0.1,2.5);
  L_Omega ~ lkj_corr_cholesky(2.);
  L_Sigma <- diag_pre_multiply(L_sigma, L_Omega);

  # Delta ~ normal(0,sigma_Delta);

  for (d in 1:D) {
    mag_int[d] ~ multi_normal_cholesky(means[d], L_Sigma);
    mag_obs[d] ~ multi_normal(mag_int[d]+gamma*k[d], mag_cov[d]);
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