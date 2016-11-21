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
  vector[5] a[6];
}

parameters {
  vector[5] c_raw;
  vector[5] alpha_raw;
  vector[5] beta_raw;
  vector<lower=0.0>[N_mags] L_sigma_raw;
  vector[5] eta_raw;

  real <lower=0> Delta_scale;
  cholesky_factor_corr[N_mags] L_Omega;

  vector[2] EW[D];
  vector[D] sivel;
  vector[N_mags] mag_int_raw[D];

  simplex[D] Delta_unit;

  vector<lower=0,upper=2.>[D] AV;
  vector<lower=1./5.5,upper=1./0.5>[D] RVinv;
  real<lower=0> AVscale;
}

transformed parameters {
  vector[5] c;
  vector[5] alpha;
  vector[5] beta;
  vector[5] eta;
  vector[N_mags] L_sigma;

  vector[D] Delta;
  vector[5] gamma;
  vector[5] rho1;
  vector[N_mags] mag_int[D];

  c = c_raw/1e2;
  alpha = alpha_raw/5e2;
  beta = beta_raw/2e2;
  eta = eta_raw/6e2;
  L_sigma = L_sigma_raw/100.;
  
  Delta = 4.*Delta_scale*(Delta_unit-1./D);


    # non-centered parameterization
  {
    matrix[5,5] L_Sigma;
    L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
    for (d in 1:D) {
      mag_int[d] = Delta[d] + c+ alpha*EW[d,1]  + beta*EW[d,2]   + eta*sivel[d] + L_Sigma * mag_int_raw[d];
    }
  }
}

model {
  real ebv;
  vector [5] AX;

  target += cauchy_lpdf(L_sigma | 0.1,0.1);
  target += lkj_corr_cholesky_lpdf(L_Omega | 4.);

  for (d in 1:D) {
    ebv  = AV[d]*RVinv[d];
    target += normal_lpdf(mag_int_raw[d]| 0, 1);
    AX = a[1]* AV[d] + a[2] * AV[d]^2
      + a[3]* ebv+ a[4] * ebv^2
      + a[5] * ebv*RVinv[d] + a[6]*(ebv*RVinv[d])^2;
    # print (AV[d], 1./RVinv[d],AX);
    target += multi_normal_lpdf(mag_obs[d] | mag_int[d] + AX, mag_cov[d]);
    target += multi_normal_lpdf(EW_obs[d] | EW[d], EW_cov[d]);
  }
  target += (normal_lpdf(sivel_obs | sivel,sivel_err));
  target += exponential_lpdf(AV | AVscale);
  target += cauchy_lpdf(AVscale | 0.1, 1);
}
