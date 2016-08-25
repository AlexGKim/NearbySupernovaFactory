data {
  int D;                // Number of supernovae
  int N_mags;
  int N_EWs;
  vector[N_mags] mag_obs[D];
  vector[N_EWs] EW_obs[D];
  matrix[N_mags, N_mags] mag_cov[D];
  matrix[N_EWs, N_EWs] EW_cov[D];
  vector[D] sivel_obs;
  vector[D] sivel_err;
}

parameters {
  vector[2] EW[D];
  vector[D] sivel;
  vector[N_mags] mag_int_raw[D];
  real <lower = 0> Delta_scale;
  real <lower = 0> k_scale;

  vector[5] c;
  vector[5] alpha;
  vector[5] beta;
  vector[5] eta;

  real gamma01;
  real gamma03;
  real gamma04;
  real gamma05;

  cholesky_factor_corr[N_mags] L_Omega;
  vector<lower=0.0, upper=0.12>[N_mags] L_sigma;

  simplex[D] Delta_unit;
  simplex[D] k_unit;
}

transformed parameters {
  vector[D] Delta;
  vector[5] gamma;
  vector[D] k;
  vector[N_mags] mag_int[D];

  gamma[1] = gamma01;
  gamma[2] = gamma03+1;
  gamma[3] = gamma03;
  gamma[4] = gamma04;
  gamma[5] = gamma05;

  Delta = Delta_scale*(Delta_unit - 1./D);
  k = k_scale*(k_unit - 1./D);

  # non-centered parameterization
  {
    matrix[5,5] L_Sigma;
    L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
    for (d in 1:D) {
      mag_int[d] = Delta[d] + c+ alpha*EW[d,1]  + beta*EW[d,2] + eta*sivel[d] + L_Sigma * mag_int_raw[d];
    }
  }
}

model {

  # Centered parameterization
  # vector[5] means[D];
  # matrix[5,5] L_Sigma;
  # for (d in 1:D) {
  #     means[d] = Delta[d] + c+ alpha*EW[d,1]  + beta*EW[d,2] + eta*sivel[d];
  # }
  # L_Sigma = diag_pre_multiply(L_sigma, L_Omega);

  target += (cauchy_lpdf(L_sigma | 0.1,0.1));
  target += (lkj_corr_cholesky_lpdf(L_Omega | 2.));

  for (d in 1:D) {
    # target += (multi_normal_cholesky_lpdf(mag_int[d] | means[d], L_Sigma));
    target += normal_lpdf(mag_int_raw[d] | 0, 1); # non-centered parameterization
    target += (multi_normal_lpdf(mag_obs[d] | mag_int[d]+gamma*k[d], mag_cov[d]));
    target += (multi_normal_lpdf(EW_obs[d] | EW[d], EW_cov[d]));
  } 
  target += (normal_lpdf(sivel_obs | sivel,sivel_err));
}