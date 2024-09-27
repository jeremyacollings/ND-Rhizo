
##### FITTING POPULATION MODELS

library(readxl)
library(tidyverse)
library(rstan)
library(loo)

# bring in raw data

fec_dat <- read_excel("NRdatMP.xlsx", sheet = "fecundity")
germ_dat <- read_excel("NRdatMP.xlsx", sheet = "germdat")
seeds <- read_excel("NRdatMP.xlsx", sheet = "seedavgerages", col_names = FALSE)

# fix case inconsistency 
fec_dat$competitor <- toupper(fec_dat$competitor)

# convert pot-level germ data to seed-level germ data

germ_dat2 <- pmap_dfr(germ_dat, 
                      function(pot_ID, nitrogen, rhizobia, 
                               focal, competitor, density, 
                               sp, baseline,
                               n.germ, n.seeded) {
                        data.frame(pot_ID = pot_ID, 
                                   nitrogen = nitrogen,
                                   rhizobia = rhizobia,
                                   focal = focal, 
                                   competitor = competitor, 
                                   density = density, 
                                   sp = sp,
                                   baseline = baseline, 
                                   germ = c( rep(1, n.germ),
                                             rep(0, n.seeded - n.germ) ) )
                      }
)

# get total seed output using mean seed counts

fec_dat$fit <- NA

fec_dat$fit <- ifelse(fec_dat$sp == "A" & fec_dat$nitrogen == "nitrogen", 
                      fec_dat$flower_count * as.numeric(seeds[1, 2]), 
                      ifelse(fec_dat$sp == "A" & fec_dat$nitrogen == "zero nitrogen", 
                             fec_dat$flower_count * as.numeric(seeds[2, 2]), fec_dat$fit ))

fec_dat$fit <- ifelse(fec_dat$sp == "B" & fec_dat$nitrogen == "nitrogen", 
                      fec_dat$flower_count * as.numeric(seeds[3, 2]), 
                      ifelse(fec_dat$sp == "B" & fec_dat$nitrogen == "zero nitrogen", 
                             fec_dat$flower_count * as.numeric(seeds[4, 2]), fec_dat$fit ))

fec_dat$fit <- ifelse(fec_dat$sp == "C" & fec_dat$nitrogen == "nitrogen", 
                      fec_dat$flower_count * as.numeric(seeds[5, 2]), 
                      ifelse(fec_dat$sp == "C" & fec_dat$nitrogen == "zero nitrogen", 
                             fec_dat$flower_count * as.numeric(seeds[6, 2]), fec_dat$fit ))

# do some prelim data viz

ggplot(data = fec_dat) + 
  geom_histogram(aes(x = fit, fill = sp), position = "identity", 
                 alpha = .6) + 
  theme_classic()

ggplot(data = fec_dat, aes(x = nitrogen, y = fit, color = rhizobia)) + 
  geom_jitter(position = position_jitterdodge()) + geom_boxplot() + facet_grid(~ sp) + 
  theme_classic()

ggplot(data = fec_dat, aes(x = rhizobia, y = fit, color = nitrogen)) + 
  geom_jitter(position = position_jitterdodge()) + geom_boxplot() + facet_grid(~ sp) + 
  theme_classic()

# these figs are simplified looks into the data and should not
# perfectly reflect the models we'll fit

fec_plots <- list()
for(foc in c("A", "B", "C")){
  fec_plots[[foc]][["A"]] <- ggplot(data = fec_dat[which(fec_dat$sp == foc & fec_dat$competitor %in% c("A", "ALONE")),], 
         aes(x = A_comp, y = fit)) + 
    geom_jitter(width = .25) + 
    stat_smooth(method = "nls",
                formula = y ~ a/(1+b*x),
                method.args = list(start = list(a = mean(fec_dat$fit[1:nrow(fec_dat)], na.rm = TRUE), b = .1)),
                se = FALSE, na.rm = TRUE) + theme_classic() +
    xlab("Competitor Density") + 
    facet_grid(nitrogen ~ rhizobia) + 
    ggtitle(paste(foc, "A", sep = "-"))
  
  fec_plots[[foc]][["B"]] <- ggplot(data = fec_dat[which(fec_dat$sp == foc & fec_dat$competitor %in% c("B", "ALONE")),], 
         aes(x = B_comp, y = fit)) + 
    geom_jitter(width = .25) + 
    stat_smooth(method = "nls",
                formula = y ~ a/(1+b*x),
                method.args = list(start = list(a = mean(fec_dat$fit[1:nrow(fec_dat)], na.rm = TRUE), b = .1)),
                se = FALSE, na.rm = TRUE) + theme_classic() +
    xlab("Competitor Density") + 
    facet_grid(nitrogen ~ rhizobia) + 
    ggtitle(paste(foc, "B", sep = "-"))
  
  fec_plots[[foc]][["C"]] <- ggplot(data = fec_dat[which(fec_dat$sp == foc & fec_dat$competitor %in% c("C", "ALONE")),], 
         aes(x = C_comp, y = fit)) + 
    geom_jitter(width = .25) + 
    stat_smooth(method = "nls",
                formula = y ~ a/(1+b*x),
                method.args = list(start = list(a = mean(fec_dat$fit[1:nrow(fec_dat)], na.rm = TRUE), b = .1)),
                se = FALSE, na.rm = TRUE) + theme_classic() +
    xlab("Competitor Density") + 
    facet_grid(nitrogen ~ rhizobia) + 
    ggtitle(paste(foc, "C", sep = "-"))
}

fec_plots

### fitting some models
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

Bev_fit4 <- stan("BevB4.stan", data=stan_dat, 
                 chains=4, cores=4)
Lot_fit <- stan("LotB.stan", data=stan_dat, 
                 chains=4, cores=4)
Lot_fit2 <- stan("LotUB.stan", data=stan_dat, 
                  chains=4, cores=4)
Rick_fit <- stan("RickB.stan", data=stan_dat, 
                 chains=4, cores=4)
Rick_fit2 <- stan("RickUB.stan", data=stan_dat, 
                 chains=4, cores=4)

loo_compare(loo(Rick_fit), loo(Rick_fit2), loo(Bev_fit), 
            loo(Lot_fit), loo(Lot_fit2))

saveRDS(Bev_fit, file = "BH_fit.rds")

# note to self... clean this up

stan_dat2 <- list(N = nrow(germ_dat2), 
                 germ = germ_dat2$germ, 
                 inds = c(length(which(germ_dat2$sp == "A" & germ_dat2$rhizobia == "sterile" & 
                                         germ_dat2$nitrogen == "zero nitrogen")), 
                          length(which(germ_dat2$sp == "A" & germ_dat2$rhizobia == "rhizobia" & 
                                         germ_dat2$nitrogen == "zero nitrogen")), 
                          length(which(germ_dat2$sp == "A" & germ_dat2$rhizobia == "sterile" & 
                                         germ_dat2$nitrogen == "nitrogen")), 
                          length(which(germ_dat2$sp == "A" & germ_dat2$rhizobia == "rhizobia" & 
                                         germ_dat2$nitrogen == "nitrogen")), 
                          length(which(germ_dat2$sp == "B" & germ_dat2$rhizobia == "sterile" & 
                                         germ_dat2$nitrogen == "zero nitrogen")), 
                          length(which(germ_dat2$sp == "B" & germ_dat2$rhizobia == "rhizobia" & 
                                         germ_dat2$nitrogen == "zero nitrogen")), 
                          length(which(germ_dat2$sp == "B" & germ_dat2$rhizobia == "sterile" & 
                                         germ_dat2$nitrogen == "nitrogen")), 
                          length(which(germ_dat2$sp == "B" & germ_dat2$rhizobia == "rhizobia" & 
                                         germ_dat2$nitrogen == "nitrogen")),
                          length(which(germ_dat2$sp == "C" & germ_dat2$rhizobia == "sterile" & 
                                         germ_dat2$nitrogen == "zero nitrogen")), 
                          length(which(germ_dat2$sp == "C" & germ_dat2$rhizobia == "rhizobia" & 
                                         germ_dat2$nitrogen == "zero nitrogen")), 
                          length(which(germ_dat2$sp == "C" & germ_dat2$rhizobia == "sterile" & 
                                         germ_dat2$nitrogen == "nitrogen")), 
                          length(which(germ_dat2$sp == "C" & germ_dat2$rhizobia == "rhizobia" & 
                                         germ_dat2$nitrogen == "nitrogen"))), 
                 A00 = which(germ_dat2$sp == "A" & germ_dat2$rhizobia == "sterile" & 
                               germ_dat2$nitrogen == "zero nitrogen"), 
                 A10 = which(germ_dat2$sp == "A" & germ_dat2$rhizobia == "rhizobia" & 
                               germ_dat2$nitrogen == "zero nitrogen"), 
                 A01 = which(germ_dat2$sp == "A" & germ_dat2$rhizobia == "sterile" & 
                               germ_dat2$nitrogen == "nitrogen"), 
                 A11 = which(germ_dat2$sp == "A" & germ_dat2$rhizobia == "rhizobia" & 
                               germ_dat2$nitrogen == "nitrogen"), 
                 B00 = which(germ_dat2$sp == "B" & germ_dat2$rhizobia == "sterile" & 
                               germ_dat2$nitrogen == "zero nitrogen"), 
                 B10 = which(germ_dat2$sp == "B" & germ_dat2$rhizobia == "rhizobia" & 
                               germ_dat2$nitrogen == "zero nitrogen"), 
                 B01 = which(germ_dat2$sp == "B" & germ_dat2$rhizobia == "sterile" & 
                               germ_dat2$nitrogen == "nitrogen"), 
                 B11 = which(germ_dat2$sp == "B" & germ_dat2$rhizobia == "rhizobia" & 
                               germ_dat2$nitrogen == "nitrogen"), 
                 C00 = which(germ_dat2$sp == "C" & germ_dat2$rhizobia == "sterile" & 
                               germ_dat2$nitrogen == "zero nitrogen"), 
                 C10 = which(germ_dat2$sp == "C" & germ_dat2$rhizobia == "rhizobia" & 
                               germ_dat2$nitrogen == "zero nitrogen"), 
                 C01 = which(germ_dat2$sp == "C" & germ_dat2$rhizobia == "sterile" & 
                               germ_dat2$nitrogen == "nitrogen"), 
                 C11 = which(germ_dat2$sp == "C" & germ_dat2$rhizobia == "rhizobia" & 
                               germ_dat2$nitrogen == "nitrogen"))

germ_ests <- stan("germ_est.stan", data = stan_dat2, chains = 4, 
                  cores = 4)

saveRDS(germ_ests, file = "germ_ests.rds")

lam_dat <- cbind.data.frame(meds = apply(as.data.frame(BH_fit)[,grepl("lam", names(as.data.frame(BH_fit)))], 
                                                               2, median), 
                            lows = apply(as.data.frame(BH_fit)[,grepl("lam", names(as.data.frame(BH_fit)))], 
                                                               2, quantile, 0.025), 
                            ups = apply(as.data.frame(BH_fit)[,grepl("lam", names(as.data.frame(BH_fit)))], 
                                                              2, quantile, 0.975), 
                            species = rep(1:3, 4), 
                            rhizo = c(rep("s", 3), rep("r", 3), rep("s", 3), rep("r", 3)), 
                            nitro = c(rep("z", 6), rep("n", 6)))

ggplot(data = lam_dat, aes(x = nitro, y = meds, color = rhizo, 
                           ymin = lows, ymax = ups)) +
  geom_point() + theme_classic(base_size = 15) + 
  facet_grid(~ species) + geom_errorbar(width = 0)

ggplot(data = lam_dat, aes(x = rhizo, y = meds, color = nitro, 
                           ymin = lows, ymax = ups)) +
  geom_point() + theme_classic(base_size = 15) + 
  facet_grid(~ species) + geom_errorbar(width = 0)

alpha_dat <- cbind.data.frame(meds = apply(as.data.frame(BH_fit)[,25:60], 2, median), 
                              lows = apply(as.data.frame(BH_fit)[,25:60], 2, quantile, 0.025), 
                              ups = apply(as.data.frame(BH_fit)[,25:60], 2, quantile, 0.975), 
                              foc = substr(names(as.data.frame(BH_fit)[,25:60]), 7, 7), 
                              comp = substr(names(as.data.frame(BH_fit)[,25:60]), 9, 9), 
                              rhizo = ifelse(substr(names(as.data.frame(BH_fit)[,25:60]), 11, 11) == 1, "s", "r"), 
                              nitro = ifelse(substr(names(as.data.frame(BH_fit)[,25:60]), 13, 13) == 1, "z", "n")) 

alpha_dat$rhizo_pair <- paste(alpha_dat$comp, alpha_dat$nitro, sep = "-")
alpha_dat$nitro_pair <- paste(alpha_dat$comp, alpha_dat$rhizo, sep = "-")

ggplot(data = alpha_dat, aes(x = rhizo, y = meds, color = nitro, 
                           ymin = lows, ymax = ups)) +
  geom_point() + 
  geom_line(aes(group = rhizo_pair)) + theme_classic(base_size = 15) + 
  facet_grid(~ foc) + geom_errorbar(width = 0)

ggplot(data = alpha_dat, aes(x = nitro, y = meds, color = rhizo, 
                             ymin = lows, ymax = ups)) +
  geom_point() + 
  geom_line(aes(group = nitro_pair)) + theme_classic(base_size = 15) + 
  facet_grid(~ foc) + geom_errorbar(width = 0)


