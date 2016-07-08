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
  vector<lower=-33, upper=-25>[N_mags] mag_int[D];

  # vector<lower=-0.4, upper=0.05>[5] c;

  real<lower=25*sqrt(5), upper=33*sqrt(5)> r_c;
  vector<lower=0, upper=pi()/2>[3] phi_c;
  real<lower=0, upper=pi()/2> phi_c_4;

  real<lower=sqrt(5)-0.1, upper=sqrt(5)+0.1> r_alpha;
  vector<lower=0, upper=pi()/2>[3] phi_alpha;
  real<lower=0, upper=pi()/2> phi_alpha_4;
  # vector<lower=-0.01, upper=0.01>[5] alpha;

  real<lower=sqrt(5)-0.1, upper=sqrt(5)+0.1> r_beta;
  # vector<lower=0, upper=pi()>[3] phi_beta;
  # real<lower=0, upper=2*pi()> phi_beta_4;
  vector<lower=0, upper=pi()/2>[3] phi_beta;
  real<lower=0, upper=pi()/2> phi_beta_4;

  cholesky_factor_corr[N_mags] L_Omega;
  vector<lower=0.0, upper = 0.25>[N_mags] L_sigma;

  vector<lower=0.55, upper=1.55>[4] gamma;
  real<lower=-2, upper=2> k[D-1];
}

transformed parameters {
  vector[5] alpha;
  vector[5] c;
  vector[5] beta;

  real dum;

  dum <- r_c;
  for (d in 1:3){
    c[d] <- dum * cos(phi_c[d]);
    dum <- dum *  sin(phi_c[d]);
  }
  c[4] <- dum * cos(phi_c_4);
  dum <- dum *  sin(phi_c_4);
  c[5] <- dum ;
  c <- -c;

  dum <- r_alpha;
  for (d in 1:3){
    alpha[d] <- dum * cos(phi_alpha[d]);
    dum <- dum *  sin(phi_alpha[d]);
  }
  alpha[4] <- dum * cos(phi_alpha_4);
  dum <- dum *  sin(phi_alpha_4);
  alpha[5] <- dum;
  alpha <- alpha-1.;

  dum <- r_beta;
  for (d in 1:3){
    beta[d] <- dum * cos(phi_beta[d]);
    dum <- dum *  sin(phi_beta[d]);
  }
  beta[4] <- dum * cos(phi_beta_4);
  dum <- dum *  sin(phi_beta_4);
  beta[5] <- dum ;
  beta <- beta-1;
}

model {
  vector[5] gamma_;
  real k_[D];
  matrix[5,5] L_Sigma;

  gamma_[1]<- gamma[1];
  gamma_[2] <- gamma[2];
  gamma_[3] <- 1.;
  gamma_[4] <- gamma[3];
  gamma_[5] <- gamma[4];

  k_[1] <- 0;
  for (d in 1:(D-1)) {
    k_[1+d] <- k[d];
  }

  L_sigma ~ cauchy(0.1,2.5);
  L_Omega ~ lkj_corr_cholesky(2.);
  L_Sigma <- diag_pre_multiply(L_sigma, L_Omega);

  for (d in 1:D) {
    mag_int[d] ~ multi_normal_cholesky(c + alpha*EW[d,1]  + beta*EW[d,2], L_Sigma);
    mag_obs[d] ~ multi_normal(mag_int[d]+gamma_*k_[d], mag_cov[d]);
    EW_obs[d] ~ multi_normal(EW[d], EW_cov[d]);
  }

  // Jacobian
  increment_log_prob(4*log(r_c));
  increment_log_prob(3*log(sin(phi_c[1])));
  increment_log_prob(2*log(sin(phi_c[2])));
  increment_log_prob(log(sin(phi_c[3])));

  increment_log_prob(4*log(r_alpha));
  increment_log_prob(3*log(sin(phi_alpha[1])));
  increment_log_prob(2*log(sin(phi_alpha[2])));
  increment_log_prob(log(sin(phi_alpha[3])));

  increment_log_prob(4*log(r_beta));
  increment_log_prob(3*log(sin(phi_beta[1])));
  increment_log_prob(2*log(sin(phi_beta[2])));
  increment_log_prob(log(sin(phi_beta[3])));
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