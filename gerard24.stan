#./gerard15 sample num_warmup=5000 num_samples=5000 data file=data.R init=init15.R output file=output15.csv refresh=1000



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

  unit_vector[5] e1;
  unit_vector[5] e2;
  unit_vector[5] e3;
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

  {
    matrix[5,5] Q;
    matrix[5,2] A;
    matrix[5,3] A2;
    matrix[5,4] A3;
    vector[5] ev1;
    vector[5] ev2;
    vector[5] ev3;
    real dp;

    for (d in 1:5){
      A[d,1] = gamma[d];
      A[d,2] = gamma1[d];
    }
    
    Q=qr_Q(A);
    Q=Q';

    ev3 = e1;
    for (d in 1:2){
      dp = dot_product(Q[d],e1);
      for(d2 in 1:5){
        ev3[d2] = ev3[d2] -  dp *Q[d,d2];
      }
    }
    ev3 = ev3/sqrt(sum(ev3 .* ev3));

    for (d in 1:5){
      A2[d,1] = gamma[d];
      A2[d,2] = gamma1[d];
      A2[d,3] = ev3[d];
    }

    Q=qr_Q(A2);
    Q=Q';
    ev1=e2;
    for (d in 1:3){
      dp = dot_product(Q[d],e2);
      for(d2 in 1:5){
        ev1[d2] = ev1[d2] -  dp *Q[d,d2];
      }
    }
    ev1 = ev1/sqrt(sum(ev1 .* ev1));

    for (d in 1:5){
      A3[d,1] = gamma[d];
      A3[d,2] = gamma1[d];
      A3[d,3] = ev3[d];
      A3[d,4] = ev1[d];
    }
    Q=qr_Q(A3);
    Q=Q';
    for(d2 in 1:5){
      ev2[d2] = Q[5,d2];
    }
    dp = dot_product(ev2, e3);
    ev2 = dp/fabs(dp) * ev2;

    # print(dot_product(gamma,ev1)," ",dot_product(gamma1,ev1));
    # print(dot_product(gamma,ev2)," ",dot_product(gamma1,ev2));
    # print(dot_product(gamma,ev3)," ",dot_product(gamma1,ev3));
    # print(dot_product(ev1,ev2)," ",dot_product(ev1,ev3)," ",dot_product(ev2,ev3));
    rho1 = rho11*ev1 + rho12*ev2 + rho13*ev3;
  }

  rho1 = -rho1*5;

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
  target += uniform_lpdf(rho1[5] | 0, 100);
  sum(R .* R) ~ cauchy(1e-3,0.5e-3);
}
