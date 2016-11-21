#./gerard11 sample num_warmup=5000 num_samples=5000 data file=data.R init=init11.R output file=output11.csv refresh=1000

data {
  int D;                // Number of supernovae
  vector[D] AV;
  vector[D] EBV;
  vector[D] RV;
  vector[5] AX[D];
}

parameters {
  # model is a1 AV + a2 AV**2 + a3 EBV + a4 EBV**2 + a5 EBV/RV + a6 (EBV/RV)**2
  vector[5]  a[9];
}

model {

  for (d in 1:D) {
    target += normal_lpdf(AX[d] | a[1]*AV[d]+ a[2] * AV[d]^2
      + a[3]*EBV[d]+ a[4] * EBV[d]^2
      + a[5] * EBV[d]/RV[d] + a[6] * (EBV[d]/RV[d])^2
      + a[7] * AV[d]*EBV[d]
      + a[8] * AV[d] * (EBV[d]/RV[d])
      + a[9] * EBV[d] * (EBV[d]/RV[d])
      , 0.01);
  }
}
