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
  vector<lower=-10, upper=50>[N_mags] mag_int[D];

  real<lower=-4.2, upper=0.1> c1;
  real<lower=-3.5, upper=0.1> c2;
  real<lower=-2.6, upper=0.1> c3;
  real<lower=-2.1, upper=0.1> c4;
  real<lower=-1.5, upper=0.1> c5;


  real<lower=0.0031-7*0.0009, upper=0.0031+8*0.0008> alpha1;
  real<lower=0.0005-7*0.0007, upper=0.0005+8*0.0007> alpha2;
  real<lower=0.0006-7*0.0006, upper=0.0006+8*0.0006> alpha3;
  real<lower=0.0007-7*0.0005, upper=0.0007+8*0.0005> alpha4;
  real<lower=0.0021-7*0.0004, upper=0.0021+8*0.0004> alpha5;


  # vector<lower=0.01, upper=0.045>[5] beta;
  real<lower=0.0345-5*0.0029, upper=0.0345+7*0.0027> beta1;
  real<lower=0.0274-5*0.0025, upper=0.0274+7*0.0022> beta2;
  real<lower=0.0274-5*0.0021, upper=0.0274+7*0.0021> beta3;
  real<lower=0.0223-5*0.0018, upper=0.0223+7*0.0018> beta4;
  real<lower=0.0213-5*0.0017, upper=0.0213+7*0.0016> beta5;

  # edges appear to give good separation, 2.3 better than 2.0
  # -2.7 too much
  # -2.3 seems too much
  # -2.2 seems too much
  real<lower=4.9882-1.5*0.3031, upper=6.9882+5*0.3399> gamma01;
  real<lower=3.0604-1.5*0.2142, upper=3.0604+5*0.2355> gamma02;
  real<lower=2.387-1.5*0.1858, upper=2.387+5*0.2009> gamma03;
  real<lower=1.7696-1.5*0.1713, upper=1.7696+5*0.1833> gamma04;

  # max of -2.2 not high enough 
  # try to get more overlap with lower gamma range.  -1.5 and -2

  real<lower=0, upper=4.9882-1.0*0.3031> rho11;
  real<lower=0, upper=3.0604-1.0*0.2142> rho12;
  real<lower=0, upper=2.387-1.0*0.1858> rho13;
  real<lower=0, upper=1.7696-1.0*0.17131> rho14;

  vector <lower=-20, upper=20.>[D] k_unit;

  # cholesky_factor_corr[N_mags] L_Omega;
  vector<lower=0.0, upper = 0.08>[N_mags] L_sigma;

  simplex[D] Delta_unit;
  real <lower = 10, upper = 35> Delta_scale;

  vector <lower=0, upper=6>[D] R;
}

transformed parameters {
  vector[D] Delta;
  vector[5] gamma;
  vector[5] alpha;
  vector[5] beta;
  vector[5] rho1;
  vector[D] k;
  vector[5] c;

  gamma[1] <- gamma01;
  gamma[2] <- 1.+gamma02;
  gamma[3] <- gamma02;
  gamma[4] <- gamma03;
  gamma[5] <- gamma04;

  alpha[1] <- alpha1;
  alpha[2] <- alpha2;
  alpha[3] <- alpha3;
  alpha[4] <- alpha4;
  alpha[5] <- alpha5;

  beta[1] <- beta1;
  beta[2] <- beta2;
  beta[3] <- beta3;
  beta[4] <- beta4;
  beta[5] <- beta5;

  rho1[1] <- rho11;
  rho1[2] <- 1+rho12;
  rho1[3] <- rho12;
  rho1[4] <- rho13;
  rho1[5] <- rho14;

  c[1] <- c1;
  c[2] <- c2;
  c[3] <- c3;
  c[4] <- c4;
  c[5] <- c5;


  Delta <- Delta_scale*(Delta_unit - 1./D);
  k <- (k_unit - mean(k_unit));
  # R <- (R_unit - mean(R_unit));
  # R <- R_unit / sqrt(sum(R .* R)/D);
}

model {
  vector[5] means[D];
  # matrix[5,5] L_Sigma;

  for (d in 1:D) {
      means[d] <- Delta[d] + c+ alpha*EW[d,1] + beta*EW[d,2] + rho1/2.*(pow(R[d],2));
  }

  increment_log_prob(cauchy_log(L_sigma, 0.,1));
  # increment_log_prob(lkj_corr_cholesky_log(L_Omega, 4.));
  # L_Sigma <- diag_pre_multiply(L_sigma, L_Omega);

  for (d in 1:D) {

    # increment_log_prob(multi_normal_cholesky_log(mag_int[d], means[d], L_Sigma));
    increment_log_prob(normal_log(mag_int[d], means[d], L_sigma));
    increment_log_prob(multi_normal_log(mag_obs[d],mag_int[d]+gamma*k[d], mag_cov[d]));
    increment_log_prob(multi_normal_log(EW_obs[d],EW[d], EW_cov[d]));

  }  
}
