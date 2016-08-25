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
  vector[N_mags] mag_int_raw[D];
  vector[5] c;
  vector[5] alpha;
  vector[5] beta;

  real <lower = 0> Delta_scale;
  real <lower = 0> k_scale;
  real k_zero;

  real gamma01;
  real gamma03;
  real gamma04;
  real gamma05;

  real gamma11;
  real gamma13;
  real gamma14;
  real gamma15;

  cholesky_factor_corr[N_mags] L_Omega;
  vector<lower=0.0, upper=0.12>[N_mags] L_sigma;
  simplex[D] k_unit;
  simplex[D] Delta_unit;

  simplex[2] prob0;
}

transformed parameters {
  vector[D] Delta;
  vector[5] gamma;
  vector[5] gamma1;
  vector[D] k;
  vector[N_mags] mag_int[D];


  gamma[1] = gamma01;
  gamma[2] = gamma03+1;
  gamma[3] = gamma03;
  gamma[4] = gamma04;
  gamma[5] = gamma05;


  gamma1[1] = gamma11;
  gamma1[2] = gamma13+1;
  gamma1[3] = gamma13;
  gamma1[4] = gamma14;
  gamma1[5] = gamma15;

  Delta = Delta_scale*(Delta_unit - 1./D);
  k = k_scale*(k_unit - k_zero);


  # non-centered parameterization
  {
    matrix[5,5] L_Sigma;
    L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
    for (d in 1:D) {
      mag_int[d] = Delta[d] + c+ alpha*EW[d,1]  + beta*EW[d,2] + L_sigma .* mag_int_raw[d];
    }
  }
}

model {
  # vector[5] means[D];
  # matrix[5,5] L_Sigma;

  real logprob0;
  real otherlog;

  # for (d in 1:D) {
  #     means[d] = Delta[d] + c+ alpha*EW[d,1]  + beta*EW[d,2];
  # }

  target += cauchy_lpdf(L_sigma | 0.08,0.1);
  # target += lkj_corr_cholesky_log(L_Omega, 2.));
  # L_Sigma = diag_pre_multiply(L_sigma, L_Omega);

  logprob0 = log(prob0[1]);
  otherlog = log(prob0[2]);
  for (d in 1:D) {
    # target += multi_normal_cholesky_log(mag_int[d], means[d], L_Sigma));
    # target += normal_lpdf(mag_int[d] | means[d], L_sigma);
    target += normal_lpdf(mag_int_raw[d] | 0, 1);

    target += log_sum_exp(logprob0+multi_normal_lpdf(mag_obs[d] | mag_int[d]+gamma*k[d], mag_cov[d]),
     otherlog+multi_normal_lpdf(mag_obs[d] | mag_int[d]+gamma1*k[d], mag_cov[d]));
    target += multi_normal_lpdf(EW_obs[d] | EW[d], EW_cov[d]);
  }  
}
