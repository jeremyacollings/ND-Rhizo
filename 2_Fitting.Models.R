
########## FITTING DEMOGRAPHIC MODELS ##########

set.seed(6)

library(tidyverse)
library(rstan)
library(parallel)
library(loo)

source("1_Data.Cleaning.R")

# Collecting data into list -----------------------------------------------

### when indexing, I'm switching the order to alphabetic
# Species A (Lupinus bicolor is now 3)
# Species B (Acmispon americanus is now 1)
# Species C (Collinsia grandiflora is now 2)

stan_dat <- list(N = nrow(fec_dat), 
                 S = length(unique(fec_dat$sp)), C = length(unique(fec_dat$sp)), 
                 R = length(unique(fec_dat$rhizobia)), 
                 r = ifelse(fec_dat$rhizobia == "sterile", 1, 2),
                 E = length(unique(fec_dat$nitrogen)), 
                 e = ifelse(fec_dat$nitrogen == "zero nitrogen", 1, 2), 
                 sp = ifelse(fec_dat$sp == "B", 1, 
                             ifelse(fec_dat$sp == "C", 2, 3)), 
                 AA = fec_dat$B_comp, CG = fec_dat$C_comp, LB = fec_dat$A_comp, 
                 fec = round(fec_dat$fit))

# Fitting candidate models ------------------------------------------------

LotB <- stan(file = file.path("Stan_Mods", "LotB.stan"), 
             data = stan_dat, chains = 2, cores = detectCores())

LotUB <- stan(file = file.path("Stan_Mods", "LotUB.stan"), 
             data = stan_dat, chains = 2, cores = detectCores())

# JC: Lotka-Volterra models aren't running...
# initializing in regions that generate negative expected values for negBinom
# I've solved this in the past with binding lambda to be > 25

BevB <- stan(file = file.path("Stan_Mods", "BevB.stan"), 
             data = stan_dat, chains = 2, cores = detectCores())

RickB <- stan(file = file.path("Stan_Mods", "RickB.stan"), 
             data = stan_dat, chains = 2, cores = detectCores())

RickUB <- stan(file = file.path("Stan_Mods", "RickUB.stan"), 
             data = stan_dat, chains = 2, cores = detectCores())

# Model comparison --------------------------------------------------------

loo_compare(#loo(LotB), loo(lotUB), 
            loo(BevB), loo(RickB), loo(RickUB))
