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

# transformed data{
#   matrix[D,D] rotation;
#   {
#     int a;
#     int b;
#     int c; 
#     a <-64072;
#     b <-46363;

#     for (i in 1:D) {
#       for (j in 1:D){
#         c <- (a+b) % 10000;
#         rotation[i,j] <- .001 * c;
#         a <- b;
#         b <- c;
#       }
#     }
#   }
#   rotation <- qr_Q(rotation);
# }

parameters {

  vector<lower=-200, upper=200>[2] EW[D];
  vector<lower=-2, upper=2>[N_mags] mag_int[D];

  # vector<lower=-0.1, upper=0.1>[5] c;
  vector<lower=-0.01, upper=0.01>[5] alpha;
  vector<lower=-0.1, upper=0.1>[5] beta;

  # real<lower=0> gamma0;
  vector<lower=0., upper=3>[4] gamma_;
  # real rho00;
  vector<lower=-.2, upper=.2>[5] rho0;
  vector<lower=-.1, upper=0.1>[5] rho1;

  # cholesky_factor_corr[N_mags] L_Omega;
  vector<lower=0.0, upper = 0.12>[N_mags] L_sigma;

  # simplex [D] Delta_simplex;
  # real <lower = 0, upper = D> Delta_scale;
  # vector [D-1] Delta_unit;
  # vector [D-1] k_unit;
  # vector [D-2] R_unit;
  # real <lower = 0, upper=D> k_scale;


  simplex[D] Delta_unit;
  real <lower = 0, upper = D/4.> Delta_scale;

  simplex[D] k_unit;
  real <lower = 0, upper = D/2.> k_scale;

  simplex[D] R_unit;

  # real<lower=0, upper=5> rho00;


  # real<lower = 0, upper = D> R_scale;

  # real <lower=-5, upper=5> R_mn;
  # vector<lower=-5,upper=5> [D] Delta;
  # real <lower=0, upper=0.2> sigma_Delta;
}

transformed parameters {
  vector[D] Delta;
  vector[5] gamma;
  vector[D] R;
  vector[D] k;


  # vector[5] rho0;

  # rho0[1] <-  rho0_[1];
  # rho0[2] <- rho0_[2];
  # rho0[3] <- 1.;
  # rho0[4] <- rho0_[3];
  # rho0[5] <- rho0_[4];

  gamma[1] <- gamma_[1];
  gamma[2] <- gamma_[2];
  gamma[3] <- 1.;
  gamma[4] <- gamma_[3];
  gamma[5] <- gamma_[4];


  Delta <- Delta_scale*(Delta_unit - 1./D);
  k <- k_scale*(k_unit - 1./D);

  R <- (R_unit - mean(R_unit));
  R <- R / sqrt(sum(R .* R)/D);

  # for (d in 1:D-1) {
  #   k[d] <- k_unit[d];
  # }
  # k[D] <- -sum(k_unit);

  # for (d in 1:D-1) {
  #   Delta[d] <- Delta_unit[d];
  # }
  # Delta[D] <- -sum(Delta_unit);

  # real m1;
  # real v1;
  # real term;

  # m1 <- sum(R_unit);
  # v1 <- sum(R_unit .* R_unit);
  # term <- 0.5 * sqrt(2*D - 2*v1 - m1*m1);

  # for (d in 1:D-2) {
  #   R[d] <- R_unit[d];
  # }
  # R[D-1] <- -0.5*m1 - term;
  # R[D] <- -0.5*m1 + term;

  # for (d in 1:D-2) {
  #   R[d] <- R_unit[d];
  # }
  # R[D-1] <- -1.;
  # R[D] <- 1.;
  # R <- rotation * R;

  # R <- R - sum(R)/D;
  # R <- R / sqrt(sum(R .* R)/D);
}

model {
  vector[5] means[D];
  # vector[5] gamma_;

  # vector[D] Delta_;
  # vector[D] R_;
  matrix[5,5] L_Sigma;

  for (d in 1:D) {
      # means[d] <- Delta[d] + c + alpha*EW[d,1]  + beta*EW[d,2] + rho0*R[d] + rho1*pow(R[d],2);
      means[d] <- Delta[d] + alpha*EW[d,1]  + beta*EW[d,2] + rho0*R[d] + rho1*pow(R[d],2);
  }

  increment_log_prob(cauchy_log(L_sigma, 0.1,2.5));
  # increment_log_prob(lkj_corr_cholesky_log(L_Omega, 2.));

  # L_Sigma <- diag_pre_multiply(L_sigma, L_Omega);


  for (d in 1:D) {

    # increment_log_prob(multi_normal_cholesky_log(mag_int[d], means[d], L_Sigma));
    increment_log_prob(normal_log(mag_int[d], means[d], L_sigma));
    increment_log_prob(multi_normal_log(mag_obs[d],mag_int[d]+gamma*k[d], mag_cov[d]));
    increment_log_prob(multi_normal_log(EW_obs[d],EW[d], EW_cov[d]));

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