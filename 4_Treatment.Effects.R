
########## Checking Treatment Effects ###########

set.seed(6)

library(tidyverse)
library(tidybayes)
library(scales)


# Bring in data -----------------------------------------------------------

output_list <- read_rds(file.path("RDS_Files", "output_list.RDS"))

# extract posterior samples from stanfit objects & put 'em in arrays
RickUB <- readRDS(file.path("RDS_Files", "fit_RickUB.RDS"))
lams_df <- as.data.frame(RickUB)[,grepl("lam", names(as.data.frame(RickUB)))]
lams_array <- array(as.numeric(unlist(lams_df)), dim = c(nrow(lams_df), 3, 2, 2))

alphas_df <- as.data.frame(RickUB)[,grepl("alpha", names(as.data.frame(RickUB)))]
alphas_array <- array(as.numeric(unlist(alphas_df)), dim = c(nrow(lams_df), 3, 3, 2, 2))

germ <- readRDS(file.path("RDS_Files", "germ.RDS"))
germ_df <- as.data.frame(germ)[,grepl("g", names(as.data.frame(germ)))]
germ_array <- array(as.numeric(unlist(germ_df)), dim = c(nrow(germ_df), 3))

mort <- readRDS(file.path("RDS_Files", "mort.RDS"))
mort_df <- as.data.frame(mort)[,grepl("m", names(as.data.frame(mort)))]
mort_array <- array(as.numeric(unlist(mort_df)), dim = c(nrow(mort_df), 3, 2, 2))

# Calculate contrasts -----------------------------------------------------

get_contrasts <- function(exp, cont, sp.dims){
  contrast <- exp - cont
  medians <- c(median(contrast, na.rm = TRUE), 
               apply(contrast, sp.dims, median, na.rm = TRUE))
  lowers <- c(quantile(contrast, na.rm = TRUE, 0.025),
              apply(contrast, sp.dims, quantile, 0.025, na.rm = TRUE))
  uppers <- c(quantile(contrast, na.rm = TRUE, 0.975),
              apply(contrast, sp.dims, quantile, 0.975, na.rm = TRUE))
  pd <- c(sum(contrast > 0, na.rm = TRUE, 0.025)/length(contrast),
          apply(contrast, sp.dims, function(x) sum(x > 0, na.rm = TRUE)/length(x)))
  cbind.data.frame(medians, lowers, uppers, pd)
}

# Population Parameters

par.df <- cbind(rbind(get_contrasts(lams_array[,,2,], lams_array[,,1,], 2), 
      get_contrasts(lams_array[,,,2], lams_array[,,,1], 2), 
      get_contrasts((lams_array[,,2,2] - lams_array[,,1,2]), 
                    (lams_array[,,2,1] - lams_array[,,1,1]), 2),
      get_contrasts(alphas_array[,,,2,], alphas_array[,,,1,], c(2,3)), 
      get_contrasts(alphas_array[,,,,2], alphas_array[,,,,1], c(2,3)), 
      get_contrasts((alphas_array[,,,2,2] - alphas_array[,,,1,2]), 
                    (alphas_array[,,,2,1] - alphas_array[,,,1,1]), c(2,3)), 
      get_contrasts(mort_array[,,2,], mort_array[,,1,], 2), 
      get_contrasts(mort_array[,,,2], mort_array[,,,1], 2), 
      get_contrasts((mort_array[,,2,2] - mort_array[,,1,2]), 
                    (mort_array[,,2,1] - mort_array[,,1,1]), 2)), 
      sp = c(rep(c("total", LETTERS[1:3]), 3), 
             rep(c("total", rep(LETTERS[1:3], 3)), 3),
             rep(c("total", LETTERS[1:3]), 3)),
      comp = c(rep("NA", 12), 
               rep(c("total", rep(LETTERS[1:3], each = 3)), 3),
               rep("NA", 12)),
      metric = rep(c("lambda", "alpha", "mort"), c(12, 30, 12)), 
      treatment = c(rep(c("rhizo", "nitro", "inter"), each = 4),
                    rep(c("rhizo", "nitro", "inter"), each = 10),
                    rep(c("rhizo", "nitro", "inter"), each = 4)))

# Godoy et al. 2014 Calculations

god.df <- cbind(rbind(get_contrasts(output_list[[1]][,,,2,], output_list[[1]][,,,1,], c(2,3)),
      get_contrasts(output_list[[1]][,,,,2], output_list[[1]][,,,,1], c(2,3)),
      get_contrasts((output_list[[1]][,,,2,2]-output_list[[1]][,,,1,2]), 
                    (output_list[[1]][,,,2,1]-output_list[[1]][,,,1,1]), c(2,3)),
      get_contrasts(output_list[[2]][,,,2,], output_list[[2]][,,,1,], c(2,3)),
      get_contrasts(output_list[[2]][,,,,2], output_list[[2]][,,,,1], c(2,3)),
      get_contrasts((output_list[[2]][,,,2,2]-output_list[[2]][,,,1,2]), 
                    (output_list[[2]][,,,2,1]-output_list[[2]][,,,1,1]), c(2,3)),
      get_contrasts(output_list[[3]][,,,2,], output_list[[3]][,,,1,], c(2,3)),
      get_contrasts(output_list[[3]][,,,,2], output_list[[3]][,,,,1], c(2,3)),
      get_contrasts((output_list[[3]][,,,2,2]-output_list[[3]][,,,1,2]), 
                    (output_list[[3]][,,,2,1]-output_list[[3]][,,,1,1]), c(2,3)),
      get_contrasts(output_list[[4]][,,,2,], output_list[[4]][,,,1,], c(2,3)),
      get_contrasts(output_list[[4]][,,,,2], output_list[[4]][,,,,1], c(2,3)),
      get_contrasts((output_list[[4]][,,,2,2]-output_list[[4]][,,,1,2]), 
                    (output_list[[4]][,,,2,1]-output_list[[4]][,,,1,1]), c(2,3))),
      sp = c(rep(c("total", rep(LETTERS[1:3], 3)), 12)),
      comp = c(rep(c("total", rep(LETTERS[1:3], each = 3)), 12)), 
      metric = rep(c("CR", "DR", "ND", "FI"), each = 30), 
      treatment = rep(rep(c("rhizo", "nitro", "inter"), each = 10), 4))


god.df <- god.df[which(god.df$sp != god.df$comp),]
god.df$pair <- apply(god.df[,c("sp", "comp")], 1, function(x) paste(sort(x), collapse = "."))
god.df <- god.df[!duplicated(paste(god.df$pair, god.df$metric, god.df$treatment)),]

# Spaak & De Laender 2021 Calculations

spaak.df <- cbind(rbind(get_contrasts(output_list[[5]][,,,2,], output_list[[1]][,,,1,], c(2,3)),
            get_contrasts(output_list[[5]][,,,,2], output_list[[1]][,,,,1], c(2,3)),
            get_contrasts((output_list[[5]][,,,2,2]-output_list[[1]][,,,1,2]), 
                          (output_list[[5]][,,,2,1]-output_list[[1]][,,,1,1]), c(2,3)),
            get_contrasts(output_list[[6]][,,,2,], output_list[[2]][,,,1,], c(2,3)),
            get_contrasts(output_list[[6]][,,,,2], output_list[[2]][,,,,1], c(2,3)),
            get_contrasts((output_list[[6]][,,,2,2]-output_list[[2]][,,,1,2]), 
                          (output_list[[6]][,,,2,1]-output_list[[2]][,,,1,1]), c(2,3))),
      sp = c(rep(c("total", rep(LETTERS[1:3], 3)), 6)),
      comp = c(rep(c("total", rep(LETTERS[1:3], each = 3)), 6)), 
      metric = rep(c("ND", "FI"), each = 30), 
      treatment = rep(rep(c("rhizo", "nitro", "inter"), each = 10), 2))

spaak.df <- spaak.df[which(spaak.df$sp != spaak.df$comp | spaak.df$sp == "total"),]

## Check significant contrasts --------------------------------------------

conf.cutoff <- .2

par.df[which(par.df$pd > 1-conf.cutoff | par.df$pd < 0+conf.cutoff),]
god.df[which(god.df$pd > 1-conf.cutoff | god.df$pd < 0+conf.cutoff),]
spaak.df[which(spaak.df$pd > 1-conf.cutoff | spaak.df$pd < 0+conf.cutoff),]

# Figures -----------------------------------------------------------------

## Parameter estimate figures ---------------------------------------------

# Intrinsic Growth Rates

RickUB %>%
  spread_draws(lam[f,r,n]) %>%
  ggplot(aes(x = as.factor(n), y = lam, color = as.factor(r))) + 
  stat_summary(fun.min = function(z) { quantile(z,0.025) },
               fun.max = function(z) { quantile(z,0.975) },
               fun = median, geom = "pointrange", 
               position = position_dodge(width = 1)) + 
  facet_wrap(~ f, labeller = labeller(f = c("1" = "AA", "2" = "CG", "3" = "LP"))) + 
  theme_classic(base_size = 15) + ylab("Intrinsic Growth Rate") + 
  scale_color_manual(name = "Rhizobia", labels = c("Unninoculated", "Innoculated"), 
                     values = c("#F5E663", "#84069D")) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Addition"))

# Species Interactions

RickUB %>%
  spread_draws(alpha[f,c,r,n]) %>%
  ggplot(aes(x = as.factor(n), y = alpha, color = as.factor(r))) + 
  stat_summary(fun.min = function(z) { quantile(z,0.025) },
               fun.max = function(z) { quantile(z,0.975) },
               fun = median, geom = "pointrange", 
               position = position_dodge(width = 1)) + 
  facet_wrap(c ~ f, labeller = labeller(c = c("1" = "cAA", "2" = "cCG", "3" = "cLP"), 
                                        f = c("1" = "fAA", "2" = "fCG", "3" = "fLP")), 
             scales = "free") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme_classic(base_size = 15) + ylab("Species Interaction") + 
  scale_color_manual(name = "Rhizobia", labels = c("Unninoculated", "Innoculated"), 
                     values = c("#F5E663", "#84069D")) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Addition"))

# Mortality Rate

mort %>%
  spread_draws(m[f,r,n]) %>%
  ggplot(aes(x = as.factor(n), y = m, color = as.factor(r))) + 
  stat_summary(fun.min = function(z) { quantile(z,0.025) },
               fun.max = function(z) { quantile(z,0.975) },
               fun = median, geom = "pointrange", 
               position = position_dodge(width = 1)) + 
  facet_wrap(~ f, labeller = labeller(f = c("1" = "AA", "2" = "CG", "3" = "LP"))) + 
  theme_classic(base_size = 15) + ylab("Mortality Rate") + 
  scale_color_manual(name = "Rhizobia", labels = c("Unninoculated", "Innoculated"), 
                     values = c("#F5E663", "#84069D")) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Addition"))

# Godoy et al. Niche Differences

cbind.data.frame(med = c(apply(output_list[["NDg"]], c(2,3,4,5), median)),
                 low = c(apply(output_list[["NDg"]], c(2,3,4,5), quantile, 0.025)),
                 up = c(apply(output_list[["NDg"]], c(2,3,4,5), quantile, 0.975)), 
                 f = rep(c("AA", "CG", "LP"), 3*2*2), 
                 c = rep(rep(c("AA", "CG", "LP"), each = 3), 2*2), 
                 r = rep(rep(c("uninoc", "inoc"), each = 3*3), 2), 
                 n = rep(c("cont", "add"), each = 3*3*2)) %>%
  filter(f != c) %>%
  ggplot(aes(x = n, y = med, ymin = low, ymax = up, color = r)) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  facet_wrap(c ~ f, labeller = labeller(c = c("AA" = "cAA", "CG" = "cCG", "LP" = "cLP"), 
                                        f = c("AA" = "fAA", "CG" = "fCG", "LP" = "fLP")), 
             scales = "free") + 
  theme_classic(base_size = 15) + ylab("Niche Differences (Godoy)") + 
  scale_color_manual(name = "Rhizobia", labels = c("Unninoculated", "Innoculated"), 
                     values = c("#F5E663", "#84069D")) + 
  scale_y_continuous(trans = pseudo_log_trans()) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Addition"))

# Godoy et al. Fitness Differences

cbind.data.frame(med = c(apply(output_list[["FIg"]], c(2,3,4,5), median)),
                 low = c(apply(output_list[["FIg"]], c(2,3,4,5), quantile, 0.025)),
                 up = c(apply(output_list[["FIg"]], c(2,3,4,5), quantile, 0.975)), 
                 f = rep(c("AA", "CG", "LP"), 3*2*2), 
                 c = rep(rep(c("AA", "CG", "LP"), each = 3), 2*2), 
                 r = rep(rep(c("uninoc", "inoc"), each = 3*3), 2), 
                 n = rep(c("cont", "add"), each = 3*3*2)) %>%
  filter(f != c) %>%
  ggplot(aes(x = n, y = med, ymin = low, ymax = up, color = r)) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  facet_wrap(c ~ f, labeller = labeller(c = c("AA" = "cAA", "CG" = "cCG", "LP" = "cLP"), 
                                        f = c("AA" = "fAA", "CG" = "fCG", "LP" = "fLP")), 
             scales = "free") + 
  theme_classic(base_size = 15) + ylab("Fitness Differences (Godoy)") + 
  scale_color_manual(name = "Rhizobia", labels = c("Unninoculated", "Innoculated"), 
                     values = c("#F5E663", "#84069D")) + 
  scale_y_continuous(trans = pseudo_log_trans()) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Addition"))

# Competitive Ratio

cbind.data.frame(med = c(apply(output_list[["CR"]], c(2,3,4,5), median)),
                 low = c(apply(output_list[["CR"]], c(2,3,4,5), quantile, 0.025)),
                 up = c(apply(output_list[["CR"]], c(2,3,4,5), quantile, 0.975)), 
                 f = rep(c("AA", "CG", "LP"), 3*2*2), 
                 c = rep(rep(c("AA", "CG", "LP"), each = 3), 2*2), 
                 r = rep(rep(c("uninoc", "inoc"), each = 3*3), 2), 
                 n = rep(c("cont", "add"), each = 3*3*2)) %>%
  filter(f != c) %>%
  ggplot(aes(x = n, y = med, ymin = low, ymax = up, color = r)) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  facet_wrap(c ~ f, labeller = labeller(c = c("AA" = "cAA", "CG" = "cCG", "LP" = "cLP"), 
                                        f = c("AA" = "fAA", "CG" = "fCG", "LP" = "fLP")), 
             scales = "free") + 
  theme_classic(base_size = 15) + ylab("Competitive Ratio") + 
  scale_color_manual(name = "Rhizobia", labels = c("Unninoculated", "Innoculated"), 
                     values = c("#F5E663", "#84069D")) + 
  scale_y_continuous(trans = pseudo_log_trans()) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Addition"))

# Demographic Ratio

cbind.data.frame(med = c(apply(output_list[["DR"]], c(2,3,4,5), median)),
                 low = c(apply(output_list[["DR"]], c(2,3,4,5), quantile, 0.025)),
                 up = c(apply(output_list[["DR"]], c(2,3,4,5), quantile, 0.975)), 
                 f = rep(c("AA", "CG", "LP"), 3*2*2), 
                 c = rep(rep(c("AA", "CG", "LP"), each = 3), 2*2), 
                 r = rep(rep(c("uninoc", "inoc"), each = 3*3), 2), 
                 n = rep(c("cont", "add"), each = 3*3*2)) %>%
  filter(f != c) %>%
  ggplot(aes(x = n, y = med, ymin = low, ymax = up, color = r)) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  facet_wrap(c ~ f, labeller = labeller(c = c("AA" = "cAA", "CG" = "cCG", "LP" = "cLP"), 
                                        f = c("AA" = "fAA", "CG" = "fCG", "LP" = "fLP")), 
             scales = "free") + 
  theme_classic(base_size = 15) + ylab("Demographic Ratio") + 
  scale_color_manual(name = "Rhizobia", labels = c("Unninoculated", "Innoculated"), 
                     values = c("#F5E663", "#84069D")) + 
  scale_y_continuous(trans = pseudo_log_trans()) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Addition"))

# Spaak & De Laender Niche Differences

cbind.data.frame(med = c(apply(output_list[["NDs"]], c(2,3,4,5), median, na.rm = TRUE)),
                 low = c(apply(output_list[["NDs"]], c(2,3,4,5), quantile, 0.025, na.rm = TRUE)),
                 up = c(apply(output_list[["NDs"]], c(2,3,4,5), quantile, 0.975, na.rm = TRUE)), 
                 f = rep(c("AA", "CG", "LP"), 3*2*2), 
                 c = rep(rep(c("AA", "CG", "LP"), each = 3), 2*2), 
                 r = rep(rep(c("uninoc", "inoc"), each = 3*3), 2), 
                 n = rep(c("cont", "add"), each = 3*3*2)) %>%
  filter(f != c) %>%
  ggplot(aes(x = n, y = med, ymin = low, ymax = up, color = r)) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  facet_wrap(c ~ f, labeller = labeller(c = c("AA" = "cAA", "CG" = "cCG", "LP" = "cLP"), 
                                        f = c("AA" = "fAA", "CG" = "fCG", "LP" = "fLP")), 
             scales = "free") + 
  theme_classic(base_size = 15) + ylab("Niche Differences (Spaak)") + 
  scale_color_manual(name = "Rhizobia", labels = c("Unninoculated", "Innoculated"), 
                     values = c("#F5E663", "#84069D")) + 
  scale_y_continuous(trans = pseudo_log_trans()) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Addition"))

# Spaak & De Laender Fitness Differences

cbind.data.frame(med = c(apply(output_list[["FIs"]], c(2,3,4,5), median, na.rm = TRUE)),
                 low = c(apply(output_list[["FIs"]], c(2,3,4,5), quantile, 0.025, na.rm = TRUE)),
                 up = c(apply(output_list[["FIs"]], c(2,3,4,5), quantile, 0.975, na.rm = TRUE)), 
                 f = rep(c("AA", "CG", "LP"), 3*2*2), 
                 c = rep(rep(c("AA", "CG", "LP"), each = 3), 2*2), 
                 r = rep(rep(c("uninoc", "inoc"), each = 3*3), 2), 
                 n = rep(c("cont", "add"), each = 3*3*2)) %>%
  filter(f != c) %>%
  ggplot(aes(x = n, y = med, ymin = low, ymax = up, color = r)) + 
  geom_pointrange(position = position_dodge(width = 1)) + 
  facet_wrap(c ~ f, labeller = labeller(c = c("AA" = "cAA", "CG" = "cCG", "LP" = "cLP"), 
                                        f = c("AA" = "fAA", "CG" = "fCG", "LP" = "fLP")), 
             scales = "free") + 
  theme_classic(base_size = 15) + ylab("Fitness Differences (Spaak)") + 
  scale_color_manual(name = "Rhizobia", labels = c("Unninoculated", "Innoculated"), 
                     values = c("#F5E663", "#84069D")) + 
  scale_y_continuous(trans = pseudo_log_trans()) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Addition"))

## Contrast figures -------------------------------------------------------

par.df %>% 
  filter(metric == "lambda") %>%
  ggplot(aes(x = sp, y = medians, ymin = lowers, ymax = uppers,
             color = treatment)) + 
  geom_pointrange(position = position_dodge(width = .5), size = 1, 
                  linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  ylab("Effect on Intrinsic Growth Rate") + 
  xlab("Species") + 
  scale_color_manual(name = "Treatment", 
                     labels = c("Interaction", "Nitrogen", "Rhizobia"), 
                     values = c("#9090A2", "#F5E663", "#84069D")) + 
  theme_classic(base_size = 15) + 
  scale_x_discrete(labels = c("AA", "CG", "LB", "Total")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "IGR.pdf"), width = 4.5, height = 4)

par.df %>% 
  filter(metric == "alpha" &
           sp != "total") %>%
  ggplot(aes(x = sp, y = medians, ymin = lowers, ymax = uppers,
             color = treatment)) + 
  geom_pointrange(position = position_dodge(width = .5), size = 1, 
                  linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  facet_wrap(~ comp, labeller = labeller(comp = c("A"="AA", 
                                                  "B"="CG",
                                                  "C"="LB"))) + 
  ylab("Effect on Interactions") + 
  xlab("Focal Species") + 
  scale_color_manual(name = "Treatment", 
                     labels = c("Interaction", "Nitrogen", "Rhizobia"), 
                     values = c("#9090A2", "#F5E663", "#84069D")) + 
  theme_classic(base_size = 15) + 
  scale_x_discrete(labels = c("AA", "CG", "LB")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "alpha.pdf"), width = 6, height = 4)

par.df %>% 
  filter(metric == "mort") %>%
  ggplot(aes(x = sp, y = medians, ymin = lowers, ymax = uppers,
             color = treatment)) + 
  geom_pointrange(position = position_dodge(width = .25), size = 1, 
                  linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  ylab("Effect on Mortality Rate") + 
  xlab("Species") + 
  scale_color_manual(name = "Treatment", 
                     labels = c("Interaction", "Nitrogen", "Rhizobia"), 
                     values = c("#9090A2", "#F5E663", "#84069D")) + 
  theme_classic(base_size = 15) + 
  scale_x_discrete(labels = c("AA", "CG", "LB", "Total")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "mort.pdf"), width = 4.5, height = 4)

god.df %>% 
  filter(metric == "ND") %>%
  ggplot(aes(x = pair, y = medians, ymin = lowers, ymax = uppers,
             color = treatment)) + 
  geom_pointrange(position = position_dodge(width = .5), size = 1, 
                  linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_y_continuous(trans = pseudo_log_trans()) + 
  ylab("Niche Differences (Godoy)") + 
  xlab("Species Pair") + 
  scale_color_manual(name = "Treatment", 
                     labels = c("Interaction", "Nitrogen", "Rhizobia"), 
                     values = c("#9090A2", "#F5E663", "#84069D")) + 
  theme_classic(base_size = 15) + 
  scale_x_discrete(labels = c("AA-CG", "AA-LB", "CG-LB")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "NDg.pdf"), width = 4.5, height = 4)

god.df %>% 
  filter(metric == "FI") %>%
  ggplot(aes(x = pair, y = medians, ymin = lowers, ymax = uppers,
             color = treatment)) + 
  geom_pointrange(position = position_dodge(width = .5), size = 1, 
                  linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_y_continuous(trans = pseudo_log_trans()) + 
  ylab("Niche Differences (Godoy)") + 
  xlab("Species Pair") + 
  scale_color_manual(name = "Treatment", 
                     labels = c("Interaction", "Nitrogen", "Rhizobia"), 
                     values = c("#9090A2", "#F5E663", "#84069D")) + 
  theme_classic(base_size = 15) + 
  scale_x_discrete(labels = c("AA-CG", "AA-LB", "CG-LB")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "FDg.pdf"), width = 4.5, height = 4)

god.df %>% 
  filter(metric == "CR") %>%
  ggplot(aes(x = pair, y = medians, ymin = lowers, ymax = uppers,
             color = treatment)) + 
  geom_pointrange(position = position_dodge(width = .5), size = 1, 
                  linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_y_continuous(trans = pseudo_log_trans()) + 
  ylab("Competitive Ratio") + 
  xlab("Species Pair") + 
  scale_color_manual(name = "Treatment", 
                     labels = c("Interaction", "Nitrogen", "Rhizobia"), 
                     values = c("#9090A2", "#F5E663", "#84069D")) + 
  theme_classic(base_size = 15) + 
  scale_x_discrete(labels = c("AA-CG", "AA-LB", "CG-LB")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "CR.pdf"), width = 4.5, height = 4)

god.df %>% 
  filter(metric == "DR") %>%
  ggplot(aes(x = pair, y = medians, ymin = lowers, ymax = uppers,
             color = treatment)) + 
  geom_pointrange(position = position_dodge(width = .5), size = 1, 
                  linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_y_continuous(trans = pseudo_log_trans()) + 
  ylab("Demographic Ratio") + 
  xlab("Species Pair") + 
  scale_color_manual(name = "Treatment", 
                     labels = c("Interaction", "Nitrogen", "Rhizobia"), 
                     values = c("#9090A2", "#F5E663", "#84069D")) + 
  theme_classic(base_size = 15) + 
  scale_x_discrete(labels = c("AA-CG", "AA-LB", "CG-LB")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "DR.pdf"), width = 4.5, height = 4)

spaak.df %>% 
  filter(metric == "ND" &
           sp != "total") %>%
  ggplot(aes(x = sp, y = medians, ymin = lowers, ymax = uppers,
             color = treatment)) + 
  geom_pointrange(position = position_dodge(width = .5), size = 1, 
                  linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  facet_wrap(~ comp, labeller = labeller(comp = c("A"="AA", 
                                                  "B"="CG",
                                                  "C"="LB"))) + 
  scale_y_continuous(trans = pseudo_log_trans()) + 
  ylab("Niche Differences (Spaak)") + 
  xlab("Focal Species") +
  scale_color_manual(name = "Treatment", 
                     labels = c("Interaction", "Nitrogen", "Rhizobia"), 
                     values = c("#9090A2", "#F5E663", "#84069D")) + 
  theme_classic(base_size = 15) + 
  scale_x_discrete(labels = c("AA", "CG", "LB")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "NDs.pdf"), width = 6, height = 4)

spaak.df %>% 
  filter(metric == "FI" &
           sp != "total") %>%
  ggplot(aes(x = sp, y = medians, ymin = lowers, ymax = uppers,
             color = treatment)) + 
  geom_pointrange(position = position_dodge(width = .5), size = 1, 
                  linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  facet_wrap(~ comp, labeller = labeller(comp = c("A"="AA", 
                                                  "B"="CG",
                                                  "C"="LB"))) + 
  scale_y_continuous(trans = pseudo_log_trans()) + 
  ylab("Fitness Differences (Spaak)") + 
  xlab("Focal Species") +
  scale_color_manual(name = "Treatment", 
                     labels = c("Interaction", "Nitrogen", "Rhizobia"), 
                     values = c("#9090A2", "#F5E663", "#84069D")) + 
  theme_classic(base_size = 15) + 
  scale_x_discrete(labels = c("AA", "CG", "LB")) + 
  theme(legend.position = "none")

ggsave(file.path("Figures", "FDs.pdf"), width = 6, height = 4)
