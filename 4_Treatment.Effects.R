
########## Checking Treatment Effects ###########

set.seed(6)

library(tidyverse)
library(scales)

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

par.df[which(par.df$pd > .8 | par.df$pd < .2),]

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

spaak.df <- spaak.df[which(spaak.df$sp != spaak.df$comp),]

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
