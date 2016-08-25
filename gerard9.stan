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
  vector[2] EW[D];
  vector[N_mags] mag_int[D];

  real <lower = 0> Delta_scale;
  real <lower = 0> k_scale;
  real k_zero;

  real c1;
  real c2;
  real c3;
  real c4;
  real c5;

  real alpha1;
  real alpha2;
  real alpha3;
  real alpha4;
  real alpha5;

  real beta1;
  real beta2;
  real beta3;
  real beta4;
  real beta5;

  real gamma01;
  real gamma03;
  real gamma04;
  real gamma05;

  real gamma11;
  real gamma12;
  real gamma13;
  real gamma14;
  real gamma15;

  # cholesky_factor_corr[N_mags] L_Omega;
  vector<lower=0.0>[N_mags] L_sigma;
  simplex[D] k_unit;
  simplex[D] Delta_unit;

  real <lower=0.0, upper=1.> prob0;
}

transformed parameters {
  vector[D] Delta;
  vector[5] gamma;
  vector[5] gamma1;
  vector[5] alpha;
  vector[5] beta;
  vector[5] c;
  vector[D] k;

  gamma[1] = gamma01;
  gamma[2] = gamma03+1;
  gamma[3] = gamma03;
  gamma[4] = gamma04;
  gamma[5] = gamma05;


  gamma1[1] = gamma11;
  gamma1[2] = gamma12;
  gamma1[3] = gamma13;
  gamma1[4] = gamma14;
  gamma1[5] = gamma15;

  alpha[1] = alpha1;
  alpha[2] = alpha2;
  alpha[3] = alpha3;
  alpha[4] = alpha4;
  alpha[5] = alpha5;

  beta[1] = beta1;
  beta[2] = beta2;
  beta[3] = beta3;
  beta[4] = beta4;
  beta[5] = beta5;

  c[1] = c1;
  c[2] = c2;
  c[3] = c3;
  c[4] = c4;
  c[5] = c5;

  Delta = Delta_scale*(Delta_unit - 1./D);
  k = k_scale*(k_unit - k_zero);
}

model {
  vector[5] means[D];
  # matrix[5,5] L_Sigma;

  real logprob0;
  real otherlog;

  for (d in 1:D) {
      means[d] = Delta[d] + c+ alpha*EW[d,1]  + beta*EW[d,2];
  }

  target += cauchy_lpdf(L_sigma | 0.08,0.2);
  # target += lkj_corr_cholesky_log(L_Omega, 2.));
  # L_Sigma = diag_pre_multiply(L_sigma, L_Omega);

  logprob0 = log(prob0);
  otherlog = log(1-prob0);
  for (d in 1:D) {
    # target += multi_normal_cholesky_log(mag_int[d], means[d], L_Sigma));
    target += normal_lpdf(mag_int[d] | means[d], L_sigma);

    target += log_sum_exp(logprob0+multi_normal_lpdf(mag_obs[d] | mag_int[d]+gamma*k[d], mag_cov[d]),
     otherlog+multi_normal_lpdf(mag_obs[d] | mag_int[d]+gamma1*k[d], mag_cov[d]));
    target += multi_normal_lpdf(EW_obs[d] | EW[d], EW_cov[d]);
  }  
}
