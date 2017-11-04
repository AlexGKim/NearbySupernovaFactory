generated quantities {

  vector[5] singular[10000];
  cholesky_factor_corr[5] L_Omega; 
  cov_matrix[5] Sigma;
  {
    matrix[5,5] L_Sigma;

    vector[5] L_sigma;
    int index;
    real dum;


    for (d in 1:10000) {
      index=1;
      while (index < 6){
        dum = cauchy_rng(0.1,0.1);
        if (dum > 0){
          L_sigma[index]=dum;
          index = index+1;
        }
      }
      // L_sigma = cauchy_rng(0.1,0.1);
      L_Omega = lkj_corr_cholesky_rng(5,.5);
      L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
      Sigma = L_Sigma * L_Sigma';
      singular[d] = singular_values(Sigma);

    }
  }

}

// ./test method=sample algorithm=fixed_param num_samples=1 num_warmup=1 output file=output1.csv
