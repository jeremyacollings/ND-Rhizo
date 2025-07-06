
########## Checking Demographic Parameters ###########

set.seed(6)

library(tidyverse)
library(tidybayes)
library(scales)
library(magrittr)

colors <- c("#FA9F42", "#2B4162", "#B02E0C")
spp <- c("AA", "CG", "LB")

source("functions.R")

# Bring in data -----------------------------------------------------------

RickUB <- readRDS(file.path("RDS_Files", "fit_RickUB.RDS"))
lams <- RickUB %>%
  spread_draws(lam[s,r,e]) %$%
  array(data = lam, dim = c(n_distinct(.draw), 
                            n_distinct(s), 
                            n_distinct(r), 
                            n_distinct(e)))

alphas <- RickUB %>%
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

mort <- readRDS(file.path("RDS_Files", "mort.RDS"))
mort_rates <- mort %>%
  spread_draws(m[s,r,e]) %$%
  array(data = m, dim = c(n_distinct(.draw), 
                            n_distinct(s), 
                            n_distinct(r), 
                            n_distinct(e)))

# Calculate contrasts -----------------------------------------------------

# Population Parameters

par.df <- cbind(rbind(get_contrasts(lams[,,2,], lams[,,1,], 2), 
                      get_contrasts(lams[,,,2], lams[,,,1], 2), 
                      get_contrasts((lams[,,2,2] - lams[,,1,2]), 
                                    (lams[,,2,1] - lams[,,1,1]), 2),
                      get_contrasts(alphas[,,,2,], alphas[,,,1,], c(2,3)), 
                      get_contrasts(alphas[,,,,2], alphas[,,,,1], c(2,3)), 
                      get_contrasts((alphas[,,,2,2] - alphas[,,,1,2]), 
                                    (alphas[,,,2,1] - alphas[,,,1,1]), c(2,3)), 
                      get_contrasts(mort_rates[,,2,], mort_rates[,,1,], 2), 
                      get_contrasts(mort_rates[,,,2], mort_rates[,,,1], 2), 
                      get_contrasts((mort_rates[,,2,2] - mort_rates[,,1,2]), 
                                    (mort_rates[,,2,1] - mort_rates[,,1,1]), 2)), 
                sp = c(rep(c("total", spp), 3), 
                       rep(c("total", rep(spp, 3)), 3),
                       rep(c("total", spp), 3)),
                comp = c(rep("NA", 12), 
                         rep(c("total", rep(spp, each = 3)), 3),
                         rep("NA", 12)),
                metric = rep(c("lambda", "alpha", "mort"), c(12, 30, 12)), 
                treatment = c(rep(c("rhizo", "nitro", "inter"), each = 4),
                              rep(c("rhizo", "nitro", "inter"), each = 10),
                              rep(c("rhizo", "nitro", "inter"), each = 4)))

conf.cutoff <- .2

par.df[which(par.df$pd > 1-conf.cutoff | par.df$pd < 0+conf.cutoff),]

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
  scale_color_manual(name = "Rhizobia", labels = c("Unninoculated", "Inoculated"), 
                     values = colors) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Fertilized"))

ggsave(file.path("Figures", "Demographic_Parameters", "lambda.pdf"), 
       width = 8, height = 6, units = "in")

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
                     values = colors) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Fertilized"))

ggsave(file.path("Figures", "Demographic_Parameters", "alphas.pdf"), 
       width = 8, height = 6, units = "in")

# Germination Rate

germ %>%
  spread_draws(g[f]) %>%
  ggplot(aes(x = as.factor(f), y = g)) + 
  stat_summary(fun.min = function(z) { quantile(z,0.025) },
               fun.max = function(z) { quantile(z,0.975) },
               fun = median, geom = "pointrange", 
               position = position_dodge(width = 1)) + 
  theme_classic(base_size = 15) + ylab("Germination Rate") +
  scale_x_discrete(name = "Species", labels = spp)

ggsave(file.path("Figures", "Demographic_Parameters", "germ.pdf"), 
       width = 8, height = 6, units = "in")

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
                     values = colors) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Control", "Fertilized"))

ggsave(file.path("Figures", "Demographic_Parameters", "mort.pdf"), 
       width = 8, height = 6, units = "in")

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
                     values = colors) + 
  theme_classic(base_size = 15) + 
  scale_x_discrete(labels = c("AA", "CG", "LB", "Total"))

ggsave(file.path("Figures", "Demographic_Parameters", "lambda_contrasts.pdf"), 
       width = 8, height = 6, units = "in")

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
                     values = colors) + 
  theme_classic(base_size = 15) + 
  scale_x_discrete(labels = c("AA", "CG", "LB"))

ggsave(file.path("Figures", "Demographic_Parameters", "alpha_contrasts.pdf"),
       width = 8, height = 6, units = "in")

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
                     values = colors) + 
  theme_classic(base_size = 15) + 
  scale_x_discrete(labels = c("AA", "CG", "LB", "Total"))

ggsave(file.path("Figures", "Demographic_Parameters", "mort_contrasts.pdf"), 
       width = 8, height = 6, units = "in")


# Network Diagrams --------------------------------------------------------

plot_network(pseudo_log(make_matrix(RickUB, 1, 1)), 
             c("AA", "CG", "LB")) + 
  ggtitle("Control, Control")

plot_network(pseudo_log(make_matrix(RickUB, 1, 2)), 
             c("AA", "CG", "LB")) + 
  ggtitle("Control, Nitrogen Addition")

plot_network(pseudo_log(make_matrix(RickUB, 2, 1)), 
             c("AA", "CG", "LB")) + 
  ggtitle("Inoculated, Control")

plot_network(pseudo_log(make_matrix(RickUB, 2, 2)),
             c("AA", "CG", "LB")) + 
  ggtitle("Inoculated, Nitrogen Addition")


# Comparing Bound and Unbound Interactions --------------------------------

RickB <- readRDS(file.path("RDS_Files", "fit_RickB.RDS"))

alpha_df <- RickUB %>%
  spread_draws(alpha[f,c,r,n])

alpha_df2 <- RickB %>%
  spread_draws(alpha[f,c,r,n])

rbind.data.frame(alpha_df, alpha_df2) %>%
  ungroup() %>%
  mutate(f = case_when(
    f == 1 ~ "AA", 
    f == 2 ~ "CG", 
    f == 3 ~ "LB"
  ), 
  c = case_when(
    c == 1 ~ "AA", 
    c == 2 ~ "CG", 
    c == 3 ~ "LB"
  )) %>%
  mutate(type = c(rep("unbound", nrow(alpha_df)), 
                  rep("bound", nrow(alpha_df2))), 
         interaction = paste(f, c, sep = "<-"), 
         r = case_when(
           r == 1 ~ "Uninoculated", 
           r == 2 ~ "Inoculated"
         ), 
         n = case_when(
           n == 1 ~ "Control", 
           n == 2 ~ "Fertilized"
         )) %>%
  group_by(r, n, type, interaction) %>%
  median_qi(alpha) %>%
  ggplot(aes(x = interaction, y = alpha, ymin = .lower, ymax = .upper, 
             color = type)) + 
  geom_pointrange(position = position_dodge(width = .25)) + 
  facet_wrap(r ~ n) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme_classic(base_size = 15) + 
  scale_color_manual(name = "Model Type", 
                     labels = c("Bound", "Unbound"), 
                     values = colors) + 
  ylab("Estimated Interaction") + 
  scale_x_discrete(name = "Interaction Pair") + 
  theme(axis.text.x = element_text(angle = 300, hjust = 0, vjust = .5))

ggsave(file.path("Figures", "Demographic_Parameters", "binding_effect.pdf"), 
       width = 8, height = 6, units = "in")
