
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

stan_dat2 <- list(N = nrow(germ_dat2), S = n_distinct(germ_dat2$sp),
                  germ = germ_dat2$germ, 
                  sp = ifelse(germ_dat2$sp == "B", 1, 
                              ifelse(germ_dat2$sp == "C", 2, 3)))

stan_dat3 <- list(N = nrow(fec_dat), S = n_distinct(fec_dat$sp), 
                  mort = ifelse(fec_dat$fit == 0, 1, 0), 
                  sp = ifelse(fec_dat$sp == "B", 1, 
                              ifelse(fec_dat$sp == "C", 2, 3)), 
                  micro = ifelse(fec_dat$rhizobia == "sterile", 1, 2), 
                  cont = ifelse(fec_dat$nitrogen == "zero nitrogen", 1, 2))

# Fitting model ------------------------------------------------

RickUB <- stan(file = file.path("Stan_Mods", "fit_RickUB.stan"), 
             data = stan_dat, chains = 4, cores = detectCores())

saveRDS(RickUB, file = file.path("RDS_Files", "fit_RickUB.RDS"))

GermMod <- stan(file = file.path("Stan_Mods", "germ.stan"), 
                data = stan_dat2, chains = 4, cores = detectCores())
 
saveRDS(GermMod, file = file.path("RDS_Files", "germ.RDS"))    

MortMod <- stan(file = file.path("Stan_Mods", "mort.stan"), 
                data = stan_dat3, chains = 4, cores = detectCores())

saveRDS(MortMod, file = file.path("RDS_Files", "mort.RDS"))  
