//
// treatment effects on germination??

data {
  int<lower=0> n; // total number of seeds
  int N[n]; // all seeds
  vector[n] spB; // is seed sp B?
  vector[n] spC; // is seed sp C?
  vector[n] Nitro; // nitrogen treatment; 0 = control, 1 = fertilized
  vector[n] Rhizo; // rhizobia treatment; 0 = sterile, 1 = innoculated
}


parameters {
  real beta_0;
  real beta_spB;
  real beta_spC;
  real beta_nitro;
  real beta_rhizo;
}

transformed parameters {
  real<lower=0,upper=1> theta[n];
  for(i in 1:n){
  theta[i] = inv_logit(beta_0 + beta_spB*spB[i] + beta_spC*spC[i] + 
  beta_nitro*Nitro[i] + beta_rhizo*Rhizo[i]);
  }
}

model {
  for(i in 1:n){
      N[i] ~ bernoulli(theta[i]);
  }
}

