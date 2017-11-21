#./gerard11 sample num_warmup=5000 num_samples=5000 data file=data.R init=init11.R output file=output11.csv refresh=1000

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
  vector[D] x1_obs;
  vector[D] x1_err;

  // matrix[5,5] rho1in_cov;
  vector[5] rho1_min;
  vector[5] rho1_max;
  matrix[5,5] rho1_ev;
}


parameters {
  vector[5] c_raw;
  vector[5] alpha_raw;
  vector[5] beta_raw;
  vector[5] eta_raw;
  vector[5] zeta;

  real<lower=0> gamma01;
  real gamma02;
  real gamma03;
  real gamma04;
  real gamma05;

  real rho11;
  real rho12;
  real rho13;
  real rho14;
  real rho15;

  // real<lower=rho1_min[1], upper=rho1_max[1]> rho11;
  // real<lower=rho1_min[2], upper=rho1_max[2]> rho12;
  // real<lower=rho1_min[3], upper=rho1_max[3]> rho13;
  // real<lower=rho1_min[4], upper=rho1_max[4]> rho14;
  // real<lower=rho1_min[5], upper=rho1_max[5]> rho15;

  real <lower=0> Delta_scale;

  // cholesky_factor_corr[N_mags] L_Omega;

  vector[2] EW[D];
  vector[D] sivel;
  vector[D] x1;

  simplex[D] Delta_unit;

  simplex[D] k_unit;
  simplex[D] R_unit;

  unit_vector[N_mags] ev;

  real<lower=0.0> ev_sig;
  vector[D] mag_int_raw;
}

transformed parameters {
  vector[5] c;
  vector[5] alpha;
  vector[5] beta;
  vector[5] eta;
  // vector[N_mags] L_sigma;

  vector[D] Delta;
  vector[D] k;
  vector[D] R;
  vector[5] gamma;
  vector[5] rho1;
  vector[N_mags] mag_int[D];

  c = c_raw/1e2;
  alpha = alpha_raw/5e2;
  beta = beta_raw/2e2;
  eta = eta_raw/6e2;
  // L_sigma = L_sigma_raw/100.;
  
  Delta = 4.*Delta_scale*(Delta_unit-1./D);
  k=(k_unit-1./D);
  R=(R_unit-1./D);

  gamma[1] = gamma01;
  gamma[2] = gamma02;
  gamma[3] = gamma03;
  gamma[4] = gamma04;
  gamma[5] = gamma05;
  gamma = gamma*5;

  // rho1[1] = rho11;
  // rho1[3] = rho13;
  // rho1[4] = rho14;
  // rho1[5] = rho15;
  // rho1[2] = rho12;
  // rho1 = rho1*5;

  for (d in 1:5){
    rho1[d]= rho11 * rho1_ev[d,1]  + rho12  * rho1_ev[d,2] + rho13 * rho1_ev[d,3] + rho14  * rho1_ev[d,4] + rho15 * rho1_ev[d,5];
  }

    # non-centered parameterization
  {
    // matrix[5,5] L_Sigma;
    // L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
    for (d in 1:D) {
      mag_int[d] = Delta[d] + c+ alpha*EW[d,1]  + beta*EW[d,2]  + zeta*x1[d]  + eta*sivel[d] + ev_sig * mag_int_raw[d] * ev;
    }
  }
}

model {
  target += cauchy_lpdf(ev_sig | 0.1,0.1);
  target += normal_lpdf(mag_int_raw| 0, 1);
  target += - D*ev_sig;
  for (d in 1:D) {
    target += multi_normal_lpdf(mag_obs[d] | mag_int[d]+gamma*k[d] + rho1*R[d], mag_cov[d]);
    target += multi_normal_lpdf(EW_obs[d] | EW[d], EW_cov[d]);
  }
  target += (normal_lpdf(sivel_obs | sivel,sivel_err));
  target += (normal_lpdf(x1_obs | x1,x1_err));
}
