// Germination rate model
data {
  int<lower=0> N; // # of germination data points
  int<lower=0> S; // # of species
  int germ[N]; // germination? yes (1) or no (0)
  int sp[N]; // sp of seed for germ data
}

parameters {
  real<lower=0,upper=1> g[S]; // germination rate per species
}

model {
  // priors
  g ~ beta(1,1);
  
  // likelihood
  for(i in 1:N){
    target += bernoulli_lpmf(germ[i] | g[sp[i]]);
  }
}

