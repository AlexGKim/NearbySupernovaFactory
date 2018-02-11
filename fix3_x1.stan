#./gerard11 sample num_warmup=5000 num_samples=5000 data file=data.R init=init11.R output file=output11.csv refresh=1000

functions{
  real covariance(vector v1, vector v2){
    return sum((v1-mean(v1)) .* (v2-mean(v2)))/num_elements(v1);
  }

  vector breakers(vector[] vs, vector delta){
    matrix[size(vs),size(vs)] covh;
    vector[size(vs)] t2;
    vector[size(vs)] ans;
    for (i in 1:size(vs)){
      covh[i,i] = covariance(vs[i],vs[i]);
      for (j in i+1:size(vs)){
        covh[i,j] = covariance(vs[i],vs[j]);
        covh[j,i] = covh[i,j];
      }
      t2[i] = covariance(vs[i],delta);
    }
    ans = inverse_spd(covh) * t2;
    // print (ans);
    // {
    //   vector[num_elements(delta)] test;
    //   test = delta;
    //   for (i in 1:size(vs)){
    //     test = test - (vs[i]-mean(vs[i]))*ans[i];
    //   }
    //   print ("begin");
    //   for (i in 1:size(vs)){
    //     print(covariance(vs[i],delta)," ",covariance(vs[i],test));
    //   }
    //   print ("end");
    // }
    return ans;
  }
}

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
  vector[5] zeta_raw;

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

  // vector[2] EW[D];
  vector[D] EW1;
  vector[D] EW2;
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
  vector[5] zeta;
  // vector[N_mags] L_sigma;

  vector[D] Delta;
  vector[D] k;
  vector[D] R;

  vector[5] gamma;
  vector[5] rho1;
  vector[5] phi;
  vector[N_mags] mag_int[D];


  c = c_raw/1e2;
  alpha = alpha_raw/5e2 ;
  beta = beta_raw/2e2;
  eta = eta_raw/6e2;
  zeta = zeta_raw;
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

  phi = ev_sig * ev ;

  {
    vector[7] dbreakers;

    real mnEW1;
    real mnEW2;
    real mnsivel;
    real mnx1;
    real mnmagintraw;
    vector[D] vs[7];

    vs[1] = EW1;
    vs[2] = EW2;
    vs[3] = sivel;
    vs[4] = x1;
    vs[5] = k;
    vs[6] = R;
    vs[7] = mag_int_raw;


    dbreakers  = breakers(vs, Delta);

    mnEW1 = mean(EW1);
    mnEW2 = mean(EW2);
    mnsivel = mean(sivel);
    mnx1 = mean(x1);
    mnmagintraw=mean(mag_int_raw);

    alpha=alpha + dbreakers[1];
    beta=beta + dbreakers[2];
    eta=eta + dbreakers[3];
    zeta=zeta + dbreakers[4];
    gamma = gamma + dbreakers[5];
    rho1 = rho1 + dbreakers[6];
    phi = phi + dbreakers[7];

    Delta  = Delta - dbreakers[1] * (EW1-mnEW1) - dbreakers[2] * (EW2-mnEW2) - dbreakers[3]*(sivel-mnsivel)
      - dbreakers[4]*(x1-mnx1) - dbreakers[5]*k - dbreakers[6]*R - dbreakers[7]*(mag_int_raw-mnmagintraw);
    c = c-dbreakers[1] * mnEW1 - dbreakers[2] * mnEW2 - dbreakers[3]*mnsivel
      - dbreakers[4]*mnx1  - dbreakers[7]*mnmagintraw;
  }
    # non-centered parameterization
  {
    // matrix[5,5] L_Sigma;
    // L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
    for (d in 1:D) {
      mag_int[d] = Delta[d] + c+ alpha*EW1[d]  + beta*EW2[d]  + zeta*x1[d]  + eta*sivel[d] + phi * mag_int_raw[d];
    }
  }
}

model {

  vector[2] tempEW;
  // real sdDelta;

  target += cauchy_lpdf(ev_sig | 0.1,0.1);
  target += normal_lpdf(mag_int_raw| 0, 1);

  for (d in 1:D) {
    tempEW[1] = EW1[d];
    tempEW[2] = EW2[d];
    target += multi_normal_lpdf(mag_obs[d] | mag_int[d]+gamma*k[d] + rho1*R[d], mag_cov[d]);
    target += multi_normal_lpdf(EW_obs[d] | tempEW, EW_cov[d]);
  }
  target += (normal_lpdf(sivel_obs | sivel,sivel_err));
  target += (normal_lpdf(x1_obs | x1,x1_err));

  // sum(Delta .* (EW1-mean(EW1))) ~ cauchy(0,1); 
  // sum(Delta .* (EW2-mean(EW2))) ~ cauchy(0,1); 
  // sum(Delta .* (sivel-mean(sivel))) ~ cauchy(0,1);   
  // sum(Delta .* (x1-mean(x1))) ~ cauchy(0,1);
  // sum(Delta .* k) ~ cauchy(0,1);
  // sum(Delta .* R) ~ cauchy(0,1);
  // sum(Delta .* (mag_int_raw-mean(mag_int_raw))) ~ cauchy(0,1);
}
