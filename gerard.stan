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
//  vector<lower=-2, upper=2>[2] EW0;
//  vector<lower=0.1, upper=1.5>[2] L_sigma_EW;
//  cholesky_factor_corr[2] L_Omega_EW;
  vector<lower=-5, upper=5>[2] EW[D];

  vector<lower=-1, upper=1>[4] c;
  vector<lower=-5, upper=5>[4] alpha;
  vector<lower=-5, upper=5>[4] beta;
  vector<lower=-5, upper=5>[3] gamma;

//  vector<lower=0.0, upper=.5>[4] L_sigma_color;
 // cholesky_factor_corr[4] L_Omega_color;
#  real<lower=-4, upper=4> k[D];
  real<lower=-2, upper=2> k[D];
#  vector<lower=0, upper=4>[4] EXY[D];
//  vector<lower=-5, upper=5>[4] colors[D];
#  vector<lower=0.01, upper = 2>[4] ebeta_inv;
  # real<lower=0.1, upper = 0.5> ebeta_inv;
}

model {
  vector[4] means[D];
  vector[4] gamma_;
  # vector[4] c_;

  //matrix[2,2] L_Sigma_EW;
  //matrix[4,4] L_Sigma_color;

  gamma_[1]<- gamma[1];
  gamma_[2] <- 1.;
  gamma_[3] <- gamma[2];
  gamma_[4] <- gamma[3];

  # c_[1]<- c[1];
  # c_[2] <- 0.;
  # c_[3] <- c[2];
  # c_[4] <- c[3];

  for (d in 1:(D-1)) {
      # means[d] <- c_ + alpha*EW[d,1]  + beta*EW[d,2] + gamma_*k[d] + EXY[d];
      means[d] <- c + alpha*EW[d,1]  + beta*EW[d,2] + gamma_*k[d];
  }

  # L_sigma_EW ~ cauchy(0,2.5);
  # L_Omega_EW ~ lkj_corr_cholesky(2.);
  # L_Sigma_EW <- diag_pre_multiply(L_sigma_EW, L_Omega_EW);
 // L_sigma_color ~ cauchy(0,0.05);
//  L_Omega_color ~ lkj_corr_cholesky(2.);
//  L_Sigma_color <- diag_pre_multiply(L_sigma_color, L_Omega_color);

  k[1] ~ normal(0,1e-8);
  for (d in 1:(D-1)) {
 //   EW[d] ~ multi_normal_cholesky(EW0,L_Sigma_EW);
  //  colors[d] ~ multi_normal_cholesky(means[d], L_Sigma_color);
    color_obs[d] ~ multi_normal(means[d], color_cov[d]);
    EW_obs[d] ~ multi_normal(EW[d], EW_cov[d]);
    # k[d] ~ exponential(1 ./ ebeta_inv);
#    EXY[d] ~ exponential(1 ./ ebeta_inv);
  }  
}

generated quantities {
  # corr_matrix[2] Omega_EW;
  //corr_matrix[4] Omega_color;

  # Omega_EW <- multiply_lower_tri_self_transpose(L_Omega_EW);
  //Omega_color <- multiply_lower_tri_self_transpose(L_Omega_color);
}
