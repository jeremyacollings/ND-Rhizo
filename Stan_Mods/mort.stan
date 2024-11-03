// Mortality rate model
data {
  int<lower=0> N; // # of mortality data points
  int<lower=0> S; // # of species
  int mort[N]; // mortality? yes (1) or no (0)
  int sp[N]; // sp of seed for mort data
  int micro[N]; // microbe treatment
  int cont[N]; // abiotic context treatment
}

parameters {
  real<lower=0,upper=1> m[S, 2, 2]; // mortality rate per species
}

model {
  // priors
  for(i in 1:S){
    for(j in 1:2){
      for(k in 1:2){
          m[i,j,k] ~ beta(1,1);
      }
    }
  }
  
  // likelihood
  for(i in 1:N){
    target += bernoulli_lpmf(mort[i] | m[sp[i], micro[i], cont[i]]);
  }
}

