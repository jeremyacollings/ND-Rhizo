
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
  real fec[N]; // fecundity
  real mean_fec[S,R,E]; // average fecundity without comp per treatment
  real sd_fec[S,R,E]; // sd  of fecundity without comp per treatment
}


parameters {
  real<lower=0> lam[S,R,E]; // pop lambda per treatment
  real<lower=0> sig[S,R,E]; // sd for lambda per treatment
  real<lower=0> alpha[S,C,R,E]; // pop alpha per comp per treatment
}

transformed parameters {
  real mu[N]; //
for(n in 1:N){
  mu[n] = lam[sp[n], r[n], e[n]]*exp(alpha[sp[n], 1, r[n], e[n]]*AA[n] + alpha[sp[n], 2, r[n], e[n]]*CG[n] + alpha[sp[n], 3, r[n], e[n]]*LB[n]);
  }
}


model {
  for(spec in 1:S){
    for(rhizo in 1:R){
    for(nitro in 1:E){
  lam[spec, rhizo, nitro] ~ normal(mean_fec[spec, rhizo, nitro], sd_fec[spec, rhizo, nitro]);
  sig[spec, rhizo, nitro] ~ uniform(.1, 100000);
    }
    }
  }
  for(fs in 1:S){
    for(cs in 1:C){
    for(rhizo in 1:R){
      for(nitro in 1:E){
      alpha[fs, cs, rhizo, nitro] ~ normal(1,1);
      }
    }
  }
  }
  for(n in 1:N){
  fec[n] ~ normal(mu[n], sig[sp[n], r[n], e[n]]);
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(fec[n] | mu[n], sig[sp[n], r[n], e[n]]);
  }
}