
library(readxl)
library(tidyverse)
library(rstan)
library(loo)
library(brms)

fec_dat <- read_excel("NRdatMP.xlsx", sheet = "fecundity")
germ_dat <- read_excel("NRdatMP.xlsx", sheet = "germdat")
seeds <- read_excel("NRdatMP.xlsx", sheet = "seedavgerages", col_names = FALSE)

fec_dat$competitor <- ifelse(fec_dat$competitor == "a", "A",
                             fec_dat$competitor)
fec_dat$competitor <- ifelse(fec_dat$competitor == "b", "B",
                             fec_dat$competitor)
fec_dat$competitor <- ifelse(fec_dat$competitor == "c", "C",
                             fec_dat$competitor)

unique(fec_dat$mortality) # seems like theres at least one NA
which(is.na(fec_dat$mortality)) # 727
fec_dat[727,] # there's a flower count... so I'm assuming this isn't a mortality
fec_dat[727,18] <- 0

# Ultimately, I suspect that 
# 1) spp have different mortality rates
# 2) nitrogen and rhizobia affect mortality
# 3) competition affects mortality
# 4) nitrogen, rhizobia and comp all interact to affect mortality
# this model is probably a bit much... so lets work up to it

# Do species have different mortality rates?

m1 <- brm(mortality ~ sp, family = "bernoulli", 
     data = fec_dat, cores = 4) # yup

# Do nitrogen or rhizobia affect mortality?

m2 <- brm(mortality ~ nitrogen + rhizobia + nitrogen*rhizobia, family = "bernoulli", 
     data = fec_dat, cores = 4) # positive interaction... ?

# but with sp?

m3 <- brm(mortality ~ sp + nitrogen + rhizobia + 
      sp*nitrogen + sp*rhizobia + sp*nitrogen*rhizobia, 
    family = "bernoulli", 
    data = fec_dat, cores = 4)

# nitrogen to vaguely decrease mortality
# Nitro:Rhizo interaction is positive overall
# C:Nitro:Rhizo is negative...

# Does competition affect mortality?

fec_dat$tot <- rowSums(fec_dat[,9:11])

# sp vary in sens to comp but not on comp effect
m4 <- brm(mortality ~ sp + tot + sp*tot, 
    family = "bernoulli", data = fec_dat, 
    cores = 4) # nope!

# sp vary in sens and comp effect
m5 <- brm(mortality ~ sp + A_comp + B_comp + C_comp + 
      sp*A_comp + sp*B_comp + sp*C_comp, 
    family = "bernoulli", data = fec_dat, 
    cores = 4)
# competing with A and C generally increases mortality
# competing with C decreases C's mortality :/

# I don't recommend running this
# for shits and gigs... 
# m6 <- brm(mortality ~ sp + A_comp + B_comp + C_comp + 
#       sp*A_comp + sp*B_comp + sp*C_comp + 
#       sp*A_comp*nitrogen + sp*B_comp*nitrogen + sp*C_comp*nitrogen + 
#       sp*A_comp*rhizobia + sp*B_comp*rhizobia + sp*C_comp*rhizobia + 
#       sp*A_comp*nitrogen*rhizobia + 
#       sp*B_comp*nitrogen*rhizobia + 
#       sp*C_comp*nitrogen*rhizobia,
#     family = "bernoulli", data = fec_dat, 
#     cores = 4)
# woot woot four-way interactions and overfitting




