//
// estimating germination rates

// The input data is a vector 'y' of length 'N'.
data {
  int N; // total number of seeds
  int germ[N]; // germ status
  int inds[12]; // indexes to keep track of... IDK how to do this in a not stupid way
  int A00[inds[1]]; // sp A, sterile, zero
  int A10[inds[2]]; // sp A, rhizo, zero
  int A01[inds[3]]; // sp A, sterile, nitro
  int A11[inds[4]]; // sp A, rhizo, nitro
  int B00[inds[5]]; // sp B, sterile, zero
  int B10[inds[6]]; // sp B, rhizo, zero
  int B01[inds[7]]; // sp B, sterile, nitro
  int B11[inds[8]]; // sp B, rhizo, nitro
  int C00[inds[9]]; // sp C, sterile, zero
  int C10[inds[10]]; // sp C, rhizo, zero
  int C01[inds[11]]; // sp C, sterile, nitro
  int C11[inds[12]]; // sp C, rhizo, nitro
}

parameters {
  real<lower=0,upper=1> rate[3, 2, 2];
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  for(i in A00){
    germ[i] ~ bernoulli(rate[1,1,1]);
  }
  
  for(i in A10){
    germ[i] ~ bernoulli(rate[1,2,1]);
  }
  
  for(i in A01){
    germ[i] ~ bernoulli(rate[1,1,2]);
  }
  
  for(i in A11){
    germ[i] ~ bernoulli(rate[1,2,2]);
  }
  
  for(i in B00){
    germ[i] ~ bernoulli(rate[2,1,1]);
  }
  
  for(i in B10){
    germ[i] ~ bernoulli(rate[2,2,1]);
  }
  
  for(i in B01){
    germ[i] ~ bernoulli(rate[2,1,2]);
  }
  
  for(i in B11){
    germ[i] ~ bernoulli(rate[2,2,2]);
  }
  
  for(i in C00){
    germ[i] ~ bernoulli(rate[3,1,1]);
  }
  
  for(i in C10){
    germ[i] ~ bernoulli(rate[3,2,1]);
  }
  
  for(i in C01){
    germ[i] ~ bernoulli(rate[3,1,2]);
  }
  
  for(i in C11){
    germ[i] ~ bernoulli(rate[3,2,2]);
  }
}

