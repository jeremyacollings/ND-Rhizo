
data {
  int N; // number of plants
  int S; // number of focal species
  int C; // number of comp species
  int R; // number of rhizobial treatments
  int r[N]; // rhizobial treatment
  int E; // number of nitrogen treatments
  int e[N]; // nitrogen treatment
  int sp[N]; // sp of focal plant
  int AA[N]; // number of CA comps
  int CG[N]; // number of CG comps
  int LB[N]; // number of PC comps
  int fec[N]; // fecundity
}


parameters {
  real<lower=0> lam[S,R,E]; // pop lambda per treatment
  real<lower=0> alpha[S,C,R,E]; // pop alpha per comp per treatment
}

transformed parameters {
  real mu[N]; //
for(n in 1:N){
  mu[n] = lam[sp[n], r[n], e[n]]*exp(alpha[sp[n], 1, r[n], e[n]]*AA[n] + alpha[sp[n], 2, r[n], e[n]]*CG[n] + alpha[sp[n], 3, r[n], e[n]]*LB[n]);
  }
}


model {
  // priors
  for(spp in 1:S){
    for(rh in 1:R){
      for(eu in 1:E){
        lam[spp,rh,eu] ~ normal(250, 100);
        for(co in 1:C){
          alpha[spp,co,rh,eu] ~ normal(0, 1);
        }
        }
        }
        }

  // likelihood
  for(n in 1:N){
  fec[n] ~ poisson(mu[n]);
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = poisson_lpmf(fec[n] | mu[n]);
  }
}