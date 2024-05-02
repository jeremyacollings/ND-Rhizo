
library(tidyverse)

fec_dat <- read_excel("NRdatMP.xlsx", sheet = "fecundity")
germ_dat <- read_excel("NRdatMP.xlsx", sheet = "germdat")
seeds <- read_excel("NRdatMP.xlsx", sheet = "seedavgerages", col_names = FALSE)

fec_dat$competitor <- ifelse(fec_dat$competitor == "a", "A",
                             fec_dat$competitor)
fec_dat$competitor <- ifelse(fec_dat$competitor == "b", "B",
                             fec_dat$competitor)
fec_dat$competitor <- ifelse(fec_dat$competitor == "c", "C",
                             fec_dat$competitor)

# alter germination data for modeling

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


# get total seed output

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


fec_dat$tot <- rowSums(fec_dat[,c("A_count", "B_count", "C_count")])

### Raw Data Figures ----- 

ggplot(data = fec_dat, aes(x = tot, y = log(fit), color = sp)) + 
  geom_point(alpha = .25) + 
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 3, b = 2)),
              se = FALSE, na.rm = TRUE) + 
  facet_wrap(nitrogen ~ rhizobia, 
             labeller = labeller(nitrogen = c("nitrogen" = "N Addition", 
                                              "zero nitrogen" = "N Control"), 
                                 rhizobia = c("rhizobia" = "Innoculated", 
                                              "sterile" = "Uninnoculated"))) + 
  xlab("Total Competitor Count") + ylab("log(Fitness)") + 
  scale_color_manual(name = "Species", labels = c("Lupinus", 
                                                  "Acmispon", 
                                                  "Collinsia"), 
                     values = c("#566E3D", "#B9A44C", "#FE9920"))

ggplot(data = germ_dat, aes(x = rhizobia, y = n.germ/n.seeded)) + 
  geom_jitter(width = .25) + 
  facet_wrap(~ sp, labeller = labeller(sp = c("A" = "Lupinus", 
                                              "B" = "Acmispon", 
                                              "C" = "Collinsia"))) + 
  geom_violin(fill = NA) + 
  xlab("Rhizobial Innoculation") + ylab("Proportion Germinated") + 
  scale_x_discrete(name = "Rhizobial Innoculation", 
                   labels = c("Innoculated", "Unninoculated")) + 
  theme_classic(base_size = 15)

### Parameter Figures -----

BH_fit <- readRDS("BH_fit.rds")
germ_ests <- readRDS("germ_ests.rds")

lams_df <- as.data.frame(BH_fit)[,1:12]
lams_array <- array(as.numeric(unlist(lams_df)), dim = c(nrow(lams_df), 3, 2, 2))

alphas_df <- as.data.frame(BH_fit)[,25:60]
alphas_array <- array(as.numeric(unlist(alphas_df)), dim = c(nrow(lams_df), 3, 3, 2, 2))

germ_df <- as.data.frame(germ_ests)[,1:12]
germ_array <- array(as.numeric(unlist(germ_df)), dim = c(nrow(germ_df), 3, 2, 2))


lam_dat <- cbind.data.frame(meds = apply(lams_df, 2, median), 
                            ups = apply(lams_df, 2, quantile, 0.975),
                            lows = apply(lams_df, 2, quantile, 0.025), 
                            rhizo = str_sub(names(lams_df), 7, 7), 
                            nitro = str_sub(names(lams_df), 9, 9), 
                            sp = str_sub(names(lams_df), 5, 5))

ggplot(data = lam_dat, aes(x = as.factor(rhizo), 
                           y = meds, ymin = lows, ymax = ups, 
                           color = as.factor(nitro), 
                           shape = as.factor(sp))) + 
  geom_line(aes(group = as.factor(paste(sp, nitro)))) + 
  geom_point(size = 3) + 
  ylab("λ") + 
  scale_x_discrete(name = "Rhizobia Treatment", 
                     labels = c("Control", "Innoculated")) + 
  theme_classic(base_size = 15) + 
  scale_color_manual(name = "Nitrogen Treatment",
                     labels = c("Control", "Addition"), 
                     values = c("#FE9920", "#566E3D")) + 
  scale_shape_manual(name = "Species", 
                     labels = c("Acmispon", "Collinsia", "Lupine"), 
                     values = c(15, 2, 16))


alpha_dat <- cbind.data.frame(meds = apply(alphas_df, 2, median), 
                              ups = apply(alphas_df, 2, quantile, 0.975),
                              lows = apply(alphas_df, 2, quantile, 0.025), 
                              rhizo = str_sub(names(alphas_df), 11, 11), 
                              nitro = str_sub(names(alphas_df), 13, 13), 
                              foc_sp = as.factor(str_sub(names(alphas_df), 7, 7)),
                              comp_sp = as.factor(str_sub(names(alphas_df), 9, 9)))

ggplot(data = alpha_dat, aes(x = as.factor(rhizo), 
                           y = meds, ymin = lows, ymax = ups, 
                           color = as.factor(nitro), 
                           shape = as.factor(comp_sp))) + 
  geom_line(aes(group = as.factor(paste(comp_sp, nitro)))) + 
  geom_point(size = 3) + 
  ylab("Sensitivity to Competition") + 
  scale_x_discrete(name = "Rhizobia Treatment", 
                   labels = c("Control", "Innoculated")) + 
  theme_classic(base_size = 15) + 
  scale_color_manual(name = "Nitrogen Treatment",
                     labels = c("Control", "Addition"), 
                     values = c("#FE9920", "#566E3D")) + 
  scale_shape_manual(name = "Competitor Species", 
                     labels = c("Acmispon", "Collinsia", "Lupine"), 
                     values = c(15, 2, 16)) + 
  facet_wrap(~ foc_sp, 
             labeller = labeller(foc_sp = c("1" = "Acmispon", 
                                   "2" = "Collinsia", 
                                   "3" = "Lupinus")))

ggplot(data = alpha_dat, aes(x = as.factor(rhizo), 
                             y = meds, ymin = lows, ymax = ups, 
                             color = as.factor(nitro), 
                             shape = as.factor(foc_sp))) + 
  geom_line(aes(group = as.factor(paste(foc_sp, nitro)))) + 
  geom_point(size = 3) + 
  ylab("Competitive Effect") + 
  scale_x_discrete(name = "Rhizobia Treatment", 
                   labels = c("Control", "Innoculated")) + 
  theme_classic(base_size = 15) + 
  scale_color_manual(name = "Nitrogen Treatment",
                     labels = c("Control", "Addition"), 
                     values = c("#FE9920", "#566E3D")) + 
  scale_shape_manual(name = "Focal Species", 
                     labels = c("Acmispon", "Collinsia", "Lupine"), 
                     values = c(15, 2, 16)) + 
  facet_wrap(~ comp_sp, 
             labeller = labeller(comp_sp = c("1" = "Acmispon", 
                                            "2" = "Collinsia", 
                                            "3" = "Lupinus")))

# Niche and Fitness Differences -----

ND <- CR1 <- CR2 <- CR <- DR1 <- DR2 <- DR <- FI <- array(NA, dim = c(nrow(lams_df), 3, 3, 2, 2))
for(i in 1:3){
  for(j in 1:3){
    ND[,i,j,,] <- 1 - sqrt((alphas_array[,i,j,,]*alphas_array[,j,i,,])/(alphas_array[,i,i,,]*alphas_array[,j,j,,]))
    CR1[,i,j,,] <- sqrt((alphas_array[,j,i,,]*alphas_array[,j,j,,])/(alphas_array[,i,i,,]*alphas_array[,i,j,,]))
    CR2[,i,j,,] <- sqrt((alphas_array[,i,i,,]*alphas_array[,i,j,,])/(alphas_array[,j,i,,]*alphas_array[,j,j,,]))
    DR1[,i,j,,] <- (germ_array[,i,,]*lams_array[,i,,])/(germ_array[,j,,]*lams_array[,j,,])
    DR2[,i,j,,] <- (germ_array[,j,,]*lams_array[,j,,])/(germ_array[,i,,]*lams_array[,i,,])
    FI[,i,j,,] <- ifelse(DR1[,i,j,,]*CR1[,i,j,,] > DR2[,i,j,,]*CR2[,i,j,,], 
                         DR1[,i,j,,]*CR1[,i,j,,], DR2[,i,j,,]*CR2[,i,j,,])
    CR[,i,j,,] <- ifelse(DR1[,i,j,,]*CR1[,i,j,,] > DR2[,i,j,,]*CR2[,i,j,,], 
                         CR1[,i,j,,], CR2[,i,j,,])
    DR[,i,j,,] <- ifelse(DR1[,i,j,,]*CR1[,i,j,,] > DR2[,i,j,,]*CR2[,i,j,,], 
                         DR1[,i,j,,], DR2[,i,j,,])
  }
}

for(i in c(1,2,3)){
  for(j in c(1,2,3)[which(1:3 > i)]){
    print(paste(i,j))
  }
  }

##### ND -----

med <- up <- low <- pair <- rhizo <- nitro <- c()
for(i in c(1,2,3)){
  for(j in c(1,2,3)[which(1:3 > i)]){
    for(k in 1:2){
      for(l in 1:2){
        med <- c(med,median(c(ND[,i,j,k,l], ND[,j,i,k,l])))
        low <- c(low,quantile(c(ND[,i,j,k,l], ND[,j,i,k,l]), 0.025))
        up <- c(up,quantile(c(ND[,i,j,k,l], ND[,j,i,k,l]), 0.975))
        pair <- c(pair, paste(i,j))
        rhizo <- c(rhizo, k)
        nitro <- c(nitro, l)
      }
    }
  }
}

ND_dat <- cbind.data.frame(med, up, low, pair, rhizo, nitro)

ggplot(data = ND_dat, aes(x = as.factor(rhizo), 
                           y = med, 
                           color = as.factor(nitro), 
                           shape = as.factor(pair))) + 
  geom_line(aes(group = as.factor(paste(pair, nitro)))) + 
  geom_point(size = 3) + 
  ylab("Niche Differences") + 
  scale_x_discrete(name = "Rhizobia Treatment", 
                   labels = c("Control", "Innoculated")) + 
  theme_classic(base_size = 15) + 
  scale_color_manual(name = "Nitrogen Treatment",
                     labels = c("Control", "Addition"), 
                     values = c("#FE9920", "#566E3D")) + 
  scale_shape_manual(name = "Species Pair", 
                     labels = c("Acmispon-Collinsia", 
                                "Acmispon-Lupine", 
                                "Collinsia-Lupine"), 
                     values = c(15, 2, 16))

ggplot(data = ND_dat, aes(x = as.factor(rhizo), 
                          y = med, ymin = low, ymax = up,
                          color = as.factor(nitro), 
                          shape = as.factor(pair))) + 
  geom_line(aes(group = as.factor(paste(pair, nitro)))) + 
  geom_point(size = 3) + 
  geom_errorbar(width = 0) +
  ylab("Niche Differences") + 
  scale_x_discrete(name = "Rhizobia Treatment", 
                   labels = c("Control", "Innoculated")) + 
  theme_classic(base_size = 15) + 
  scale_color_manual(name = "Nitrogen Treatment",
                     labels = c("Control", "Addition"), 
                     values = c("#FE9920", "#566E3D")) + 
  scale_shape_manual(name = "Species Pair", 
                     labels = c("Acmispon-Collinsia", 
                                "Acmispon-Lupine", 
                                "Collinsia-Lupine"), 
                     values = c(15, 2, 16))

ggplot(data = ND_dat, aes(x = as.factor(rhizo), 
                          y = low, 
                          color = as.factor(nitro), 
                          shape = as.factor(pair))) + 
  geom_line(aes(group = as.factor(paste(pair, nitro)))) + 
  geom_point(size = 3) + 
  ylab("Niche Differences") + 
  scale_x_discrete(name = "Rhizobia Treatment", 
                   labels = c("Control", "Innoculated")) + 
  theme_classic(base_size = 15) + 
  scale_color_manual(name = "Nitrogen Treatment",
                     labels = c("Control", "Addition"), 
                     values = c("#FE9920", "#566E3D")) + 
  scale_shape_manual(name = "Species Pair", 
                     labels = c("Acmispon-Collinsia", 
                                "Acmispon-Lupine", 
                                "Collinsia-Lupine"), 
                     values = c(15, 2, 16))

ggplot(data = ND_dat, aes(x = as.factor(rhizo), 
                          y = up, 
                          color = as.factor(nitro), 
                          shape = as.factor(pair))) + 
  geom_line(aes(group = as.factor(paste(pair, nitro)))) + 
  geom_point(size = 3) + 
  ylab("Niche Differences") + 
  scale_x_discrete(name = "Rhizobia Treatment", 
                   labels = c("Control", "Innoculated")) + 
  theme_classic(base_size = 15) + 
  scale_color_manual(name = "Nitrogen Treatment",
                     labels = c("Control", "Addition"), 
                     values = c("#FE9920", "#566E3D")) + 
  scale_shape_manual(name = "Species Pair", 
                     labels = c("Acmispon-Collinsia", 
                                "Acmispon-Lupine", 
                                "Collinsia-Lupine"), 
                     values = c(15, 2, 16))

##### FI -----

med <- up <- low <- pair <- rhizo <- nitro <- c()
for(i in c(1,2,3)){
  for(j in c(1,2,3)[which(1:3 > i)]){
    for(k in 1:2){
      for(l in 1:2){
        med <- c(med,median(c(FI[,i,j,k,l], FI[,j,i,k,l])))
        low <- c(low,quantile(c(FI[,i,j,k,l], FI[,j,i,k,l]), 0.025))
        up <- c(up,quantile(c(FI[,i,j,k,l], FI[,j,i,k,l]), 0.975))
        pair <- c(pair, paste(i,j))
        rhizo <- c(rhizo, k)
        nitro <- c(nitro, l)
      }
    }
  }
}

FI_dat <- cbind.data.frame(med, up, low, pair, rhizo, nitro)

ggplot(data = FI_dat, aes(x = as.factor(rhizo), 
                          y = med, 
                          color = as.factor(nitro), 
                          shape = as.factor(pair))) + 
  geom_line(aes(group = as.factor(paste(pair, nitro)))) + 
  geom_point(size = 3) + 
  ylab("Fitness Inequalities") + 
  scale_x_discrete(name = "Rhizobia Treatment", 
                   labels = c("Control", "Innoculated")) + 
  theme_classic(base_size = 15) + 
  scale_color_manual(name = "Nitrogen Treatment",
                     labels = c("Control", "Addition"), 
                     values = c("#FE9920", "#566E3D")) + 
  scale_shape_manual(name = "Species Pair", 
                     labels = c("Acmispon-Collinsia", 
                                "Acmispon-Lupine", 
                                "Collinsia-Lupine"), 
                     values = c(15, 2, 16))

ggplot(data = FI_dat, aes(x = as.factor(rhizo), 
                          y = med, ymin = low, ymax = up,
                          color = as.factor(nitro), 
                          shape = as.factor(pair))) + 
  geom_line(aes(group = as.factor(paste(pair, nitro)))) + 
  geom_point(size = 3) + 
  geom_errorbar(width = 0) +
  ylab("Fitness Inequalities") + 
  scale_x_discrete(name = "Rhizobia Treatment", 
                   labels = c("Control", "Innoculated")) + 
  theme_classic(base_size = 15) + 
  scale_color_manual(name = "Nitrogen Treatment",
                     labels = c("Control", "Addition"), 
                     values = c("#FE9920", "#566E3D")) + 
  scale_shape_manual(name = "Species Pair", 
                     labels = c("Acmispon-Collinsia", 
                                "Acmispon-Lupine", 
                                "Collinsia-Lupine"), 
                     values = c(15, 2, 16))

ggplot(data = FI_dat, aes(x = as.factor(rhizo), 
                          y = low, 
                          color = as.factor(nitro), 
                          shape = as.factor(pair))) + 
  geom_line(aes(group = as.factor(paste(pair, nitro)))) + 
  geom_point(size = 3) + 
  ylab("Fitness Inequalities") + 
  scale_x_discrete(name = "Rhizobia Treatment", 
                   labels = c("Control", "Innoculated")) + 
  theme_classic(base_size = 15) + 
  scale_color_manual(name = "Nitrogen Treatment",
                     labels = c("Control", "Addition"), 
                     values = c("#FE9920", "#566E3D")) + 
  scale_shape_manual(name = "Species Pair", 
                     labels = c("Acmispon-Collinsia", 
                                "Acmispon-Lupine", 
                                "Collinsia-Lupine"), 
                     values = c(15, 2, 16))

ggplot(data = FI_dat, aes(x = as.factor(rhizo), 
                          y = up, 
                          color = as.factor(nitro), 
                          shape = as.factor(pair))) + 
  geom_line(aes(group = as.factor(paste(pair, nitro)))) + 
  geom_point(size = 3) + 
  ylab("Fitness Inequalities") + 
  scale_x_discrete(name = "Rhizobia Treatment", 
                   labels = c("Control", "Innoculated")) + 
  theme_classic(base_size = 15) + 
  scale_color_manual(name = "Nitrogen Treatment",
                     labels = c("Control", "Addition"), 
                     values = c("#FE9920", "#566E3D")) + 
  scale_shape_manual(name = "Species Pair", 
                     labels = c("Acmispon-Collinsia", 
                                "Acmispon-Lupine", 
                                "Collinsia-Lupine"), 
                     values = c(15, 2, 16))

### Coexistence Results -----

pair <- prob <- med_ND <- med_FI <- rhizo <- nitro <- c()

for(i in c(1,2,3)){
  for(j in c(1,2,3)[which(1:3 > i)]){
    for(k in 1:2){
      for(l in 1:2){
        prob <- c(prob, (sum(ND[,i,j,k,l] > FI[,i,j,k,l] - 1) + 
                    sum(ND[,j,i,k,l] > FI[,j,i,k,l] - 1))/
          length(c(ND[,i,j,k,l], ND[,j,i,k,l])))
        med_ND <- c(med_ND,median(c(ND[,i,j,k,l], ND[,j,i,k,l])))
        med_FI <- c(med_FI,median(c(FI[,i,j,k,l], FI[,j,i,k,l])))
        pair <- c(pair, paste(i,j))
        rhizo <- c(rhizo, k)
        nitro <- c(nitro, l)
      }
    }
  }
}

NDs <- FIs <- c()
for(i in c(1,2,3)){
  for(j in c(1,2,3)){
    for(k in 1:2){
      for(l in 1:2){
        NDs <- c(NDs, ND[,i,j,k,l], ND[,j,i,k,l])
        FIs <- c(FIs, FI[,i,j,k,l], FI[,j,i,k,l])
      }
    }
  }
}

sum(NDs > FIs)

cx_dat <- cbind.data.frame(pair, prob, med_ND, med_FI, rhizo, nitro)

find_hull <- function(df) df[chull(df$med_ND, df$med_FI), ]
hulls <- plyr::ddply(cx_dat, "pair", find_hull)

ggplot(data = cx_dat, aes(x = med_ND, y = med_FI, 
                          group = paste(pair, rhizo, nitro), 
                          shape = paste(rhizo, nitro), 
                          color = pair)) + 
  geom_point(size = 3) + 
  geom_polygon(data = hulls, aes(group = pair),
               alpha = 0.5, fill = NA, size = 1) + 
  stat_function(fun = function(x) -x+1, color = "black") + 
  stat_function(fun = function(x) 1/(-x+1), color = "black") + 
  #ylim(.97,8) + 
  scale_shape_manual(name = "Treatment", 
                     labels = c("Unninoculated, N Control", 
                                "Unninoculated, N Addition", 
                                "Inoculated, N Control", 
                                "Inoculated, N Addition"), 
                     values = c(1, 16, 2, 17)) + 
  scale_color_manual(name = "Species Pair", 
                     labels = c("Acmispon-Collinsia", 
                                "Acmispon-Lupine", 
                                "Collinsia-Lupine"), 
                     values = c("#566E3D", "#FE9920", "#B9A44C")) + 
  xlab("Niche Differences") + ylab("Fitness Inequalities") + 
  theme_classic(base_size = 15)


ggplot(data = cx_dat, aes(x = as.factor(rhizo), 
                          y = prob, 
                          color = as.factor(nitro), 
                          shape = as.factor(pair))) + 
  geom_line(aes(group = as.factor(paste(pair, nitro)))) + 
  geom_point(size = 3) + 
  ylab("Coexistence Probability") + 
  scale_x_discrete(name = "Rhizobia Treatment", 
                   labels = c("Control", "Innoculated")) + 
  theme_classic(base_size = 15) + 
  scale_color_manual(name = "Nitrogen Treatment",
                     labels = c("Control", "Addition"), 
                     values = c("#FE9920", "#566E3D")) + 
  scale_shape_manual(name = "Species Pair", 
                     labels = c("Acmispon-Collinsia", 
                                "Acmispon-Lupine", 
                                "Collinsia-Lupine"), 
                     values = c(15, 2, 16))




