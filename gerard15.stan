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
  vector[5] c_raw;
  vector[5] alpha_raw;
  vector[5] beta_raw;
  vector<lower=0.0>[N_mags] L_sigma_raw;
  vector[5] eta_raw;

  real<lower=0> gamma01;
  real gamma02;
  real gamma03;
  real gamma04;
  real gamma05;

  real<upper=0> gamma11;
  real gamma12;
  real gamma13;
  real gamma14;
  real gamma15;

  real rho11;
  real rho12;
  real rho13;
  # real rho14;
  # real rho15;

  real <lower=0> Delta_scale;
  cholesky_factor_corr[N_mags] L_Omega;

  vector[2] EW[D];
  vector[D] sivel;
  vector[N_mags] mag_int_raw[D];

  simplex[D] Delta_unit;

  simplex[D] k_unit;
  simplex[D] k1_unit;

  simplex[D] R_unit;
}

transformed parameters {
  vector[5] c;
  vector[5] alpha;
  vector[5] beta;
  vector[5] eta;
  vector[N_mags] L_sigma;

  vector[D] Delta;
  vector[D] k;
  vector[D] k1;
  vector[D] R;
  vector[5] gamma;
  vector[5] gamma1;
  vector[5] rho1;
  vector[N_mags] mag_int[D];

  c = c_raw/1e2;
  alpha = alpha_raw/5e2;
  beta = beta_raw/2e2;
  eta = eta_raw/6e2;
  L_sigma = L_sigma_raw/100.;
  
  Delta = 4.*Delta_scale*(Delta_unit-1./D);
  k=(k_unit-1./D);
  k1=(k1_unit-1./D);
  R=(R_unit-1./D);

  gamma[1] = gamma01;
  gamma[2] = gamma02;
  gamma[3] = gamma03;
  gamma[4] = gamma04;
  gamma[5] = gamma05;
  gamma = gamma*5;

  gamma1[1] = gamma11;
  gamma1[2] = gamma12;
  gamma1[3] = gamma13;
  gamma1[4] = gamma14;
  gamma1[5] = gamma15;
  gamma1 = gamma1*5;

  rho1[1] = (-rho11/2+rho12/2+rho13);
  rho1[2] = (rho11/2-rho12-rho13/2);
  rho1[3] = (rho11/2+rho12-rho13/2);
  rho1[4] = (rho11/2+rho12/2-rho13);
  rho1[5] = (-rho11+rho12/2+rho13/2);
  rho1 = rho1*5;
  {
    vector[2] y;
    matrix[2,2] A;
    vector[2] x;
    y[1] = dot_product(rho1, gamma);
    y[2] = dot_product(rho1, gamma1);
    A[1,1] = dot_product(gamma,gamma);
    A[2,2] = dot_product(gamma1,gamma1);
    A[1,2] = dot_product(gamma,gamma1);
    A[2,1] = A[1,2];
    x = inverse(A) * y;
    rho1= rho1 - x[1]*gamma - x[2]*gamma1;
  }

  # {
  #   matrix[5,5] Q;
  #   matrix[5,2] A;
  #   for (d in 1:5){
  #     A[d,1] = gamma[d];
  #     A[d,2] = gamma1[d];
  #   }
  #   Q = qr_Q(A);
  #   for (d in 1:5){
  #     rho1[d] = rho11*Q[d,3]+rho12*Q[d,4]+rho13*Q[d,5];
  #   }
  #   print (rho11,rho12,rho13);
  #   rho1=rho1*5;
  # }

    # non-centered parameterization
  {
    matrix[5,5] L_Sigma;
    L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
    for (d in 1:D) {
      mag_int[d] = Delta[d] + c+ alpha*EW[d,1]  + beta*EW[d,2]  + rho1*R[d]  + eta*sivel[d] + L_Sigma * mag_int_raw[d];
    }
  }
}

model {
  target += cauchy_lpdf(L_sigma | 0.1,0.1);
  target += lkj_corr_cholesky_lpdf(L_Omega | 4.);

  for (d in 1:D) {
    target += normal_lpdf(mag_int_raw[d]| 0, 1);
    target += multi_normal_lpdf(mag_obs[d] | mag_int[d]+gamma*k[d]+gamma1*k1[d], mag_cov[d]);
    target += multi_normal_lpdf(EW_obs[d] | EW[d], EW_cov[d]);
  }
  target += (normal_lpdf(sivel_obs | sivel,sivel_err));
}
