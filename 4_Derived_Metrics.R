
########## Checking Demographic Parameters ###########

set.seed(6)

library(tidyverse)
library(tidybayes)
library(scales)
library(magrittr)
library(reshape2)

colors <- c("#FA9F42", "#2B4162", "#B02E0C")
spp <- c("AA", "CG", "LB")

source("functions.R")

# Bring in data -----------------------------------------------------------

RickB <- readRDS(file.path("RDS_Files", "fit_RickB.RDS"))
lams <- RickB %>%
  spread_draws(lam[s,r,e]) %$%
  array(data = lam, dim = c(n_distinct(.draw), 
                            n_distinct(s), 
                            n_distinct(r), 
                            n_distinct(e)))

alphas <- RickB %>%
  spread_draws(alpha[s,c, r,e]) %$%
  array(data = alpha, dim = c(n_distinct(.draw), 
                              n_distinct(s), 
                              n_distinct(c),
                              n_distinct(r), 
                              n_distinct(e)))

germ <- readRDS(file.path("RDS_Files", "germ.RDS"))
germ_rates <- germ %>%
  spread_draws(g[s]) %$%
  array(data = g, dim = c(n_distinct(.draw), 
                          n_distinct(s)))

# Calculate Niche and Fitness Differences ---------------------------------

# initialize arrays to store metrics 
ND <- CR1 <- CR2 <- CR <- DR1 <- DR2 <- DR <- FI <- array(NA, dim = dim(alphas))
CX <- winners <- array(NA, dim = dim(alphas))
for(i in 1:3){
  for(j in 1:3){
    CR1[,i,j,,] <- sqrt((alphas[,j,i,,]*alphas[,j,j,,])/(alphas[,i,i,,]*alphas[,i,j,,]))
    CR2[,i,j,,] <- sqrt((alphas[,i,i,,]*alphas[,i,j,,])/(alphas[,j,i,,]*alphas[,j,j,,]))
    DR1[,i,j,,] <- (germ_rates[,i]*lams[,i,,])/(germ_rates[,j]*lams[,j,,])
    DR2[,i,j,,] <- (germ_rates[,j]*lams[,j,,])/(germ_rates[,i]*lams[,i,,])
    FI[,i,j,,] <- ifelse(DR1[,i,j,,]*CR1[,i,j,,] > DR2[,i,j,,]*CR2[,i,j,,], 
                         DR1[,i,j,,]*CR1[,i,j,,], DR2[,i,j,,]*CR2[,i,j,,])
    CR[,i,j,,] <- ifelse(DR1[,i,j,,]*CR1[,i,j,,] > DR2[,i,j,,]*CR2[,i,j,,],
                         CR1[,i,j,,], CR2[,i,j,,])
    DR[,i,j,,] <- ifelse(DR1[,i,j,,]*CR1[,i,j,,] > DR2[,i,j,,]*CR2[,i,j,,],
                         DR1[,i,j,,], DR2[,i,j,,])
    ND[,i,j,,] <- 1 - sqrt((alphas[,i,j,,]*alphas[,j,i,,])/
                             (alphas[,i,i,,]*alphas[,j,j,,]))
    
    # need to write code to say 1. is ND = 1 & FI = 0... then NA
    # 2. do ND and FI suggest CX... then 1
    # 3. do ND and FI suggest EX... then 2
    # 4. do ND and FI suggest PE... then 3
    
    # then a seperate array of the winners... 
    # if first array != 2... then NA
    # if first array == 2...
    # who won? 1, 2, or 3? 
    
    CX[,i,j,,] <- ifelse(ND[,i,j,,] == 1 & FI[,i,j,,] == 0, NA, 
                         ifelse(1 - ND[,i,j,,] < FI[,i,j,,] &
                                  FI[,i,j,,] < 1 / (1 - ND[,i,j,,]), 1, # coexistence
                                ifelse(1 - ND[,i,j,,] < FI[,i,j,,] &
                                         FI[,i,j,,] > 1 / (1 - ND[,i,j,,]), 2, # exclusions
                                       ifelse(1 - ND[,i,j,,] > FI[,i,j,,] &
                                                FI[,i,j,,] > 1 / (1 - ND[,i,j,,]), 3, NA)))) # priority effects
    
    winners[,i,j,,] <- ifelse(CX[,i,j,,] %in% c(NA, 1, 3), NA, 
                              ifelse(CX[,i,j,,] == 2 &
                                       DR1[,i,j,,]*CR1[,i,j,,] > DR2[,i,j,,]*CR2[,i,j,,], i, 
                                     ifelse(CX[,i,j,,] == 2 &
                                              DR1[,i,j,,]*CR1[,i,j,,] < DR2[,i,j,,]*CR2[,i,j,,], j, NA)))
  }
}

saveRDS(list("CR" = CR, "DR" = DR, "ND" = ND, "FI" = FI, 
             "CX" = CX, "winners" = winners), 
        file.path("RDS_Files", "derived_quantities.rds"))

# Calculate Contrasts -----------------------------------------------------

conts <- cbind(rbind(get_contrasts(CR[,,,2,], CR[,,,1,], c(2,3)),
                      get_contrasts(CR[,,,,2], CR[,,,,1], c(2,3)),
                      get_contrasts((CR[,,,2,2]-CR[,,,1,2]), 
                                    (CR[,,,2,1]-CR[,,,1,1]), c(2,3)),
                      get_contrasts(DR[,,,2,], DR[,,,1,], c(2,3)),
                      get_contrasts(DR[,,,,2], DR[,,,,1], c(2,3)),
                      get_contrasts((DR[,,,2,2]-DR[,,,1,2]), 
                                    (DR[,,,2,1]-DR[,,,1,1]), c(2,3)),
                      get_contrasts(ND[,,,2,], ND[,,,1,], c(2,3)),
                      get_contrasts(ND[,,,,2], ND[,,,,1], c(2,3)),
                      get_contrasts((ND[,,,2,2]-ND[,,,1,2]), 
                                    (ND[,,,2,1]-ND[,,,1,1]), c(2,3)),
                      get_contrasts(FI[,,,2,], FI[,,,1,], c(2,3)),
                      get_contrasts(FI[,,,,2], FI[,,,,1], c(2,3)),
                      get_contrasts((FI[,,,2,2]-FI[,,,1,2]), 
                                    (FI[,,,2,1]-FI[,,,1,1]), c(2,3))),
                sp = c(rep(c("total", rep(spp, 3)), 12)),
                comp = c(rep(c("total", rep(spp, each = 3)), 12)), 
                metric = rep(c("CR", "DR", "ND", "FI"), each = 30), 
                treatment = rep(rep(c("rhizo", "nitro", "inter"), each = 10), 4))


conts <- conts[which(conts$sp != conts$comp),]
conts$pair <- apply(conts[,c("sp", "comp")], 1, function(x) paste(sort(x), collapse = "."))
conts <- conts[!duplicated(paste(conts$pair, conts$metric, conts$treatment)),]

conf.cutoff <- .2

conts[which(conts$pd > 1-conf.cutoff | conts$pd < 0+conf.cutoff),]

# Figures -----------------------------------------------------------------

## Parameter estimate figures ---------------------------------------------

CR_df <- array2df(CR, "CR")
DR_df <- array2df(DR, "DR")
ND_df <- array2df(ND, "ND")
FI_df <- array2df(FI, "FI")

# Competitive Ratio

CR_df %>%
  ggplot(aes(x = as.factor(nitro), y = CR, color = as.factor(rhizo))) + 
  stat_summary(fun.min = function(z) { quantile(z,0.025) },
               fun.max = function(z) { quantile(z,0.975) },
               fun = median, geom = "pointrange", 
               position = position_dodge(width = 1)) + 
  facet_wrap(~ pair, labeller = labeller(pair = c("1.2" = "AA.CG", "1.3" = "AA.LP", "2.3" = "CG.LP"))) + 
  theme_classic(base_size = 15) + ylab("Competitive Ratio") + 
  scale_color_manual(name = "Rhizobia", labels = c("Unninoculated", "Inoculated"), 
                     values = colors) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Fertilized")) + 
  theme(axis.text.x = element_text(angle = 300, hjust = 0, vjust = .5))

ggsave(file.path("Figures", "Coexistence_Metrics", "CR.pdf"), 
       width = 8, height = 6, units = "in")

# Demographic Ratio

DR_df %>%
  ggplot(aes(x = as.factor(nitro), y = DR, color = as.factor(rhizo))) + 
  stat_summary(fun.min = function(z) { quantile(z,0.025) },
               fun.max = function(z) { quantile(z,0.975) },
               fun = median, geom = "pointrange", 
               position = position_dodge(width = 1)) + 
  facet_wrap(~ pair, labeller = labeller(pair = c("1.2" = "AA.CG", "1.3" = "AA.LP", "2.3" = "CG.LP"))) + 
  theme_classic(base_size = 15) + ylab("Demographic Ratio") + 
  scale_color_manual(name = "Rhizobia", labels = c("Unninoculated", "Inoculated"), 
                     values = colors) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Fertilized")) + 
  theme(axis.text.x = element_text(angle = 300, hjust = 0, vjust = .5))

ggsave(file.path("Figures", "Coexistence_Metrics", "DR.pdf"), 
       width = 8, height = 6, units = "in")

# Niche Differences

ND_df %>%
  ggplot(aes(x = as.factor(nitro), y = ND, color = as.factor(rhizo))) + 
  stat_summary(fun.min = function(z) { quantile(z,0.025) },
               fun.max = function(z) { quantile(z,0.975) },
               fun = median, geom = "pointrange", 
               position = position_dodge(width = 1)) + 
  facet_wrap(~ pair, labeller = labeller(pair = c("1.2" = "AA.CG", "1.3" = "AA.LP", "2.3" = "CG.LP"))) + 
  theme_classic(base_size = 15) + ylab("Niche Differences") + 
  scale_color_manual(name = "Rhizobia", labels = c("Unninoculated", "Inoculated"), 
                     values = colors) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Fertilized")) + 
  theme(axis.text.x = element_text(angle = 300, hjust = 0, vjust = .5))

ggsave(file.path("Figures", "Coexistence_Metrics", "ND.pdf"), 
       width = 8, height = 6, units = "in")

# Fitness Inequalities

FI_df %>%
  ggplot(aes(x = as.factor(nitro), y = FI, color = as.factor(rhizo))) + 
  stat_summary(fun.min = function(z) { quantile(z,0.025) },
               fun.max = function(z) { quantile(z,0.975) },
               fun = median, geom = "pointrange", 
               position = position_dodge(width = 1)) + 
  facet_wrap(~ pair, labeller = labeller(pair = c("1.2" = "AA.CG", "1.3" = "AA.LP", "2.3" = "CG.LP"))) + 
  theme_classic(base_size = 15) + ylab("Fitness Inequalities") + 
  scale_color_manual(name = "Rhizobia", labels = c("Unninoculated", "Inoculated"), 
                     values = colors) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Fertilized")) + 
  theme(axis.text.x = element_text(angle = 300, hjust = 0, vjust = .5))

ggsave(file.path("Figures", "Coexistence_Metrics", "FI.pdf"), 
       width = 8, height = 6, units = "in")

## Contrast figures -------------------------------------------------------

# Competitive Ratio

conts %>% 
  filter(metric == "CR") %>%
  ggplot(aes(x = pair, y = medians, ymin = lowers, ymax = uppers,
             color = treatment)) + 
  geom_pointrange(position = position_dodge(width = .5), size = 1, 
                  linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  ylab("Effect on Competitive Ratio") + 
  xlab("Species Pair") + 
  scale_color_manual(name = "Treatment", 
                     labels = c("Interaction", "Nitrogen", "Rhizobia"), 
                     values = colors) + 
  theme_classic(base_size = 15) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "Coexistence_Metrics", "CR_Contrasts.pdf"), width = 4.5, height = 4)

# Demographic Ratio

conts %>% 
  filter(metric == "DR") %>%
  ggplot(aes(x = pair, y = medians, ymin = lowers, ymax = uppers,
             color = treatment)) + 
  geom_pointrange(position = position_dodge(width = .5), size = 1, 
                  linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  ylab("Effect on Demographic Ratio") + 
  xlab("Species Pair") + 
  scale_color_manual(name = "Treatment", 
                     labels = c("Interaction", "Nitrogen", "Rhizobia"), 
                     values = colors) + 
  theme_classic(base_size = 15) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "Coexistence_Metrics", "DR_Contrasts.pdf"), width = 4.5, height = 4)

# Niche Differences

conts %>% 
  filter(metric == "ND") %>%
  ggplot(aes(x = pair, y = medians, ymin = lowers, ymax = uppers,
             color = treatment)) + 
  geom_pointrange(position = position_dodge(width = .5), size = 1, 
                  linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  ylab("Effect on Niche Differences") + 
  xlab("Species Pair") + 
  scale_color_manual(name = "Treatment", 
                     labels = c("Interaction", "Nitrogen", "Rhizobia"), 
                     values = colors) + 
  theme_classic(base_size = 15) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "Coexistence_Metrics", "ND_Contrasts.pdf"), width = 4.5, height = 4)

# Fitness Inequalities

conts %>% 
  filter(metric == "FI") %>%
  ggplot(aes(x = pair, y = medians, ymin = lowers, ymax = uppers,
             color = treatment)) + 
  geom_pointrange(position = position_dodge(width = .5), size = 1, 
                  linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  ylab("Effect on Fitness Inequalities") + 
  xlab("Species Pair") + 
  scale_color_manual(name = "Treatment", 
                     labels = c("Interaction", "Nitrogen", "Rhizobia"), 
                     values = colors) + 
  theme_classic(base_size = 15) +  
  theme(legend.position = "none")

ggsave(file.path("Figures", "Coexistence_Metrics", "FI_Contrasts.pdf"), width = 4.5, height = 4)
