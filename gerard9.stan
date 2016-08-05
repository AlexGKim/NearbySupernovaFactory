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
  vector<lower=-2.5, upper=2.5>[N_mags] mag_int[D];

  vector<lower=-0.5, upper=0.15>[5] c;
  # vector<lower=-0.002, upper=0.006>[5] alpha;
  real<lower=-0.002, upper=0.003> alpha1;
  real<lower=-0.004, upper=-0.0007> alpha2;
  real<lower=-0.0027, upper=-0.0008> alpha3;
  real<lower=-0.003, upper=-0.001> alpha4;
  real<lower=0.0005, upper=0.0045> alpha5;

  # vector<lower=0.01, upper=0.045>[5] beta;
  real<lower=0.006, upper=0.021> beta1;
  real<lower=0.000, upper=0.013> beta2;
  real<lower=0.0025, upper=0.0085> beta3;
  real<lower=-0.002, upper=0.005> beta4;
  real<lower=0.015, upper=0.03> beta5;

  # vector<lower=1., upper=6>[4] gamma_;
  real<lower=1.6, upper=2.2> gamma01;
  real<lower=1.8, upper=4.3> gamma02;
  real<lower=-0.9, upper=-.5> gamma03;
  real<lower=-2, upper=-.8> gamma04;

  vector<lower=0.1, upper=4>[4] gamma1_;

  cholesky_factor_corr[N_mags] L_Omega;
  vector<lower=0.0, upper = 0.15>[N_mags] L_sigma;

  simplex[D] Delta_unit;
  real <lower = 0, upper = 3*D/4.> Delta_scale;

  real <lower=0.0, upper=1.> prob0;

  vector <lower=-0.5, upper=1.>[D] k;
}

transformed parameters {
  vector[D] Delta;
  vector[5] gamma;
  vector[5] gamma1;
  vector[5] alpha;
  vector[5] beta;
  # vector[D] k;

  gamma[1] <- gamma01+gamma02;
  gamma[2] <- 1.+gamma02;
  gamma[3] <- gamma02;
  gamma[4] <- gamma03+gamma02;
  gamma[5] <- gamma04+gamma02;

  gamma1[1] <- gamma1_[1];
  gamma1[2] <- 1.+gamma1_[2];
  gamma1[3] <- gamma1_[2];
  gamma1[4] <- gamma1_[3];
  gamma1[5] <- gamma1_[4];

  alpha[1] <- alpha1+alpha5;
  alpha[2] <- alpha2+alpha5;
  alpha[3] <- alpha3+alpha5;
  alpha[4] <- alpha4+alpha5;
  alpha[5] <- alpha5;

  beta[1] <- beta1+beta5;
  beta[2] <- beta2+beta5;
  beta[3] <- beta3+beta5;
  beta[4] <- beta4+beta5;
  beta[5] <- beta5;

  Delta <- Delta_scale*(Delta_unit - 1./D);
  # k <- k_unit - mean(k_unit);
}

model {
  vector[5] means[D];
  matrix[5,5] L_Sigma;

  real logprob0;
  real otherlog;

  for (d in 1:D) {
      # means[d] <- Delta[d] + c + alpha*EW[d,1]  + beta*EW[d,2];
      means[d] <- Delta[d] + c+ alpha*EW[d,1]  + beta*EW[d,2];
  }

  increment_log_prob(cauchy_log(L_sigma, 0.1,2.5));
  increment_log_prob(lkj_corr_cholesky_log(L_Omega, 2.));
  L_Sigma <- diag_pre_multiply(L_sigma, L_Omega);

  logprob0 <- log(prob0);
  otherlog <- log(1-prob0);
  for (d in 1:D) {
    increment_log_prob(multi_normal_cholesky_log(mag_int[d], means[d], L_Sigma));
    # increment_log_prob(normal_log(mag_int[d], means[d], L_sigma));

    increment_log_prob(log_sum_exp(logprob0+multi_normal_log(mag_obs[d],mag_int[d]+gamma*k[d], mag_cov[d]),
     otherlog+multi_normal_log(mag_obs[d],mag_int[d]+gamma1*k[d], mag_cov[d])));
    increment_log_prob(multi_normal_log(EW_obs[d],EW[d], EW_cov[d]));


  }  
}