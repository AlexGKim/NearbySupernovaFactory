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
  vector[5] c_raw;
  vector[5] alpha_raw;
  vector[5] beta_raw;
  vector<lower=0.0>[N_mags] L_sigma_raw;

  real<lower=0> gamma01;
  real<lower=0>  gamma05;
  real<lower=0>  gamma02;
  real<lower=0>  gamma03;
  real<lower=0>  gamma04;

  real<lower=0>  rho11;
  real rho12;
  real rho13;
  real rho14;
  real rho15;

  real <lower=0> Delta_scale;
  # cholesky_factor_corr[N_mags] L_Omega;

  vector[2] EW[D];
  vector[N_mags] mag_int_raw[D];

  simplex[D] Delta_unit;

  simplex[D] k_unit;
  simplex[D] R_unit;


  # real <lower=0> R_scale;
  # real <lower=0> k_scale;

}

transformed parameters {
  vector[5] c;
  vector[5] alpha;
  vector[5] beta;
  vector[N_mags] L_sigma;

  vector[D] Delta;
  vector[D] k;
  vector[D] R;
  vector[5] gamma;
  vector[5] rho1;
  vector[N_mags] mag_int[D];

  c = c_raw/1e2;
  alpha = alpha_raw/5e2;
  beta = beta_raw/2e2;
  L_sigma = L_sigma_raw/100.;

  Delta = 10*Delta_scale*(Delta_unit-1./D);
  # k=k_scale*(k_unit-1./D);
  # R=R_scale*(R_unit-1./D);

  k=(k_unit-1./D);
  R=(R_unit-1./D);

  gamma[1] = gamma01;
  gamma[3] = gamma03;
  gamma[4] = gamma04;
  gamma[5] = gamma05;
  gamma[2] = gamma02;
  gamma = gamma*5;

  rho1[1] = rho11;
  rho1[3] = rho13;
  rho1[4] = rho14;
  rho1[5] = rho15;
  rho1[2] = rho12;
  rho1 = rho1*5;

    # non-centered parameterization
  {
    # matrix[5,5] L_Sigma;
    # L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
    for (d in 1:D) {
      mag_int[d] = Delta[d] + c+ alpha*EW[d,1]  + beta*EW[d,2]  + rho1*R[d] + L_sigma .* mag_int_raw[d];
    }
  }
}

model {
  # vector[5] means[D];
  # matrix[5,5] L_Sigma;

  # for (d in 1:D) {
  #     means[d] = Delta[d] + c+ alpha*EW[d,1] + beta*EW[d,2] + rho1*R[d];
  #           # means[d] = Delta[d] + c+ alpha*EW[d,1] + beta*EW[d,2] + rho1/2.*(pow(R[d],2));
  # }
  # L_Sigma = diag_pre_multiply(L_sigma, L_Omega);

  target += cauchy_lpdf(L_sigma | 0.1,0.1);
  # target += lkj_corr_cholesky_lpdf(L_Omega | 4.);

  for (d in 1:D) {
    target += normal_lpdf(mag_int_raw[d]| 0, 1);
    target += multi_normal_lpdf(mag_obs[d] | mag_int[d]+gamma*k[d], mag_cov[d]);
    target += multi_normal_lpdf(EW_obs[d] | EW[d], EW_cov[d]);
  }  
}
