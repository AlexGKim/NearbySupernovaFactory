#limit beta to 3 sigma from 5

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
  vector<lower=-1.9, upper=2.1>[N_mags] mag_int[D];

  vector<lower=-0.1, upper=0.16>[5] c;
  # vector<lower=-0.002, upper=0.006>[5] alpha;

  real<lower=0, upper=0.0031+5*0.0008> alpha1;
  real<lower=0.00, upper=0.0005+5*0.0007> alpha2;
  real<lower=-0.002, upper=0.0006+5*0.0006> alpha3;
  real<lower=-0.002, upper=0.0007+5*0.0005> alpha4;
  real<lower=-0.002, upper=0.0021+5*0.0004> alpha5;

  # vector<lower=0.01, upper=0.045>[5] beta;
  real<lower=0.0345-5*0.0029, upper=0.0345+5*0.0027> beta1;
  real<lower=0.0274-5*0.0025, upper=0.0274+5*0.0022> beta2;
  real<lower=0.0274-5*0.0021, upper=0.0274+5*0.0021> beta3;
  real<lower=0.0223-5*0.0018, upper=0.0223+5*0.0018> beta4;
  real<lower=0.0213-5*0.0017, upper=0.0213+5*0.0016> beta5;

  # vector<lower=1., upper=6>[4] gamma_;
  real<lower=4.9882-3*0.3031, upper=4.9882+4*0.3399> gamma01;
  real<lower=3.0604-3*0.2142, upper=3.0604+4*0.2355> gamma02;
  real<lower=2.387-3*0.1858, upper=2.387+4*0.2009> gamma03;
  real<lower=1.7696-3*0.1713, upper=1.7696+4*0.1833> gamma04;

  real<lower=1.2, upper=3.> gamma11;
  real<lower=-0.1, upper=1.4> gamma12;
  real<lower=-0.5, upper=1.> gamma13;
  real<lower=-1, upper=.6> gamma14;

  # cholesky_factor_corr[N_mags] L_Omega;
  vector<lower=0.0, upper = 0.08>[N_mags] L_sigma;

  simplex[D] Delta_unit;
  real <lower = 15, upper = 45> Delta_scale;

  real <lower=0.0, upper=1.> prob0;

  vector <lower=-0.4, upper=.8>[D] k;
}

transformed parameters {
  vector[D] Delta;
  vector[5] gamma;
  vector[5] gamma1;
  vector[5] alpha;
  vector[5] beta;
  # vector[D] k;

  gamma[1] <- gamma01;
  gamma[2] <- 1.+gamma02;
  gamma[3] <- gamma02;
  gamma[4] <- gamma03;
  gamma[5] <- gamma04;

  gamma1[1] <- gamma11;
  gamma1[2] <- 1.+gamma12;
  gamma1[3] <- gamma12;
  gamma1[4] <- gamma13;
  gamma1[5] <- gamma14;

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

  Delta <- Delta_scale*(Delta_unit - 1./D);
  # k <- k_unit - mean(k_unit);
}

model {
  vector[5] means[D];
  # matrix[5,5] L_Sigma;

  real logprob0;
  real otherlog;

  for (d in 1:D) {
      # means[d] <- Delta[d] + c + alpha*EW[d,1]  + beta*EW[d,2];
      means[d] <- Delta[d] + c+ alpha*EW[d,1]  + beta*EW[d,2];
  }

  increment_log_prob(cauchy_log(L_sigma, 0.1,2.5));
  # increment_log_prob(lkj_corr_cholesky_log(L_Omega, 2.));
  # L_Sigma <- diag_pre_multiply(L_sigma, L_Omega);

  logprob0 <- log(prob0);
  otherlog <- log(1-prob0);
  for (d in 1:D) {
    # increment_log_prob(multi_normal_cholesky_log(mag_int[d], means[d], L_Sigma));
    increment_log_prob(normal_log(mag_int[d], means[d], L_sigma));

    increment_log_prob(log_sum_exp(logprob0+multi_normal_log(mag_obs[d],mag_int[d]+gamma*k[d], mag_cov[d]),
     otherlog+multi_normal_log(mag_obs[d],mag_int[d]+gamma1*k[d], mag_cov[d])));
    increment_log_prob(multi_normal_log(EW_obs[d],EW[d], EW_cov[d]));


  }  
}