
data {
  int N; // number of plants
  int S; // number of focal species
  int C; // number of comp species
  int R; // number of rhizobial treatments
  int r[N]; // rhizobial treatment
  int E; // number of nitrogen treatments
  int eu[N]; // nitrogen treatment
  int sp[N]; // sp of focal plant
  int AA[N]; // number of CA comps
  int CG[N]; // number of CG comps
  int LB[N]; // number of PC comps
  int fec[N]; // fecundity
}


parameters {
  real<lower=0> lam[S,R,E]; // pop lambda per treatment
  real<lower=0> phi[S]; // dispersion for lambda per species
  real free_alpha[S,C,R,E]; // pop alpha per comp per treatment
  real<lower=0> intra_alpha[S,R,E]; // intraspecific comp
}

transformed parameters {
  real alpha[S,C,R,E]; // pop alpha per comp per treatment
  real mu[N]; // expected fitness
  
  // make alpha array
  for(spp in 1:S){
    for(com in 1:C){
      for(rhi in 1:R){
        for(eut in 1:E){
          if(spp == com){
            alpha[spp, com, rhi, eut] = intra_alpha[spp, rhi, eut];
          }
          else{
            alpha[spp, com, rhi, eut] = free_alpha[spp, com, rhi, eut];
          }
        }
      }
    }
  }
  
  // calculate expected fitness
for(n in 1:N){
  mu[n] = lam[sp[n], r[n], eu[n]]*exp(-(alpha[sp[n], 1, r[n], eu[n]]*AA[n] + 
  alpha[sp[n], 2, r[n], eu[n]]*CG[n] + 
  alpha[sp[n], 3, r[n], eu[n]]*LB[n]));
  }
}


model {
  // priors
  for(spp in 1:S){
    for(rh in 1:R){
      for(eut in 1:E){
        lam[spp,rh,eut] ~ normal(250, 100);
        intra_alpha[spp,rh,eut] ~ normal(0, 1);
        for(co in 1:C){
          free_alpha[spp,co,rh,eut] ~ normal(0, 1);
        }
        }
        }
        }

phi ~ exponential(1); 
  
  // likelihood
  for(n in 1:N){
  fec[n] ~ neg_binomial_2(mu[n], phi[sp[n]]);
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = neg_binomial_2_lpmf(fec[n] | mu[n], phi[sp[n]]);
  }
}