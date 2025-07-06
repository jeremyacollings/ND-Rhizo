
########## EXPLORATORY DATA VISUALIZATION ##########

set.seed(6)

colors <- c("#FA9F42", "#2B4162", "#B02E0C")

library(tidyverse)

source("1_Data.Cleaning.R")

fec_dat$sp2 <- ifelse(fec_dat$sp == "A", "LB", 
                        ifelse(fec_dat$sp == "B", "AA", "CG"))


# Fitness Plots -----------------------------------------------------------

## Fitness ~ Treatments ---------------------------------------------------

# all

fec_dat %>%
  ggplot(aes(x = nitrogen, y = log(fit), color = rhizobia)) + 
  geom_jitter(position = position_jitterdodge()) + geom_boxplot() + facet_grid(~ sp2) + 
  theme_classic() + 
  scale_color_manual(values = colors, name = "Rhizobia", 
                     labels = c("Inoculated", "Uninoculated")) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Added", "Control")) + 
  ylab("log(Fitness)")

ggsave(file.path("Figures", "Raw", "Fitness_All1.pdf"), 
       width = 8, height = 6, units = "in")

fec_dat %>%
  ggplot(aes(x = rhizobia, y = log(fit), color = nitrogen)) + 
  geom_jitter(position = position_jitterdodge()) + geom_boxplot() + facet_grid(~ sp2) + 
  theme_classic() + 
  scale_color_manual(values = colors, name = "Nitrogen", 
                     labels = c("Fertilized", "Control")) + 
  scale_x_discrete(name = "Rhizobia", labels = c("Inoculated", "Uninoculated")) + 
  ylab("log(Fitness)")

ggsave(file.path("Figures", "Raw", "Fitness_All2.pdf"), 
       width = 8, height = 6, units = "in")

# just alones

fec_dat %>%
  filter(density == "alone") %>%
  ggplot(aes(x = nitrogen, y = log(fit), color = rhizobia)) + 
  geom_jitter(position = position_jitterdodge()) + geom_boxplot() + facet_grid(~ sp2) + 
  theme_classic() + 
  scale_color_manual(values = colors, name = "Rhizobia", 
                     labels = c("Inoculated", "Uninoculated")) + 
  scale_x_discrete(name = "Nitrogen", labels = c("Fertilized", "Control")) + 
  ylab("log(Fitness)")

ggsave(file.path("Figures", "Raw", "Fitness_Alones1.pdf"), 
       width = 8, height = 6, units = "in")

fec_dat %>%
  filter(density == "alone") %>%
  ggplot(aes(x = rhizobia, y = log(fit), color = nitrogen)) + 
  geom_jitter(position = position_jitterdodge()) + geom_boxplot() + facet_grid(~ sp2) + 
  theme_classic() + 
  scale_color_manual(values = colors, name = "Nitrogen", 
                     labels = c("Fertilized", "Control")) + 
  scale_x_discrete(name = "Rhizobia", labels = c("Inoculated", "Uninoculated")) + 
  ylab("log(Fitness)")

ggsave(file.path("Figures", "Raw", "Fitness_Alones2.pdf"), 
       width = 8, height = 6, units = "in")


## Fitness ~ Competition --------------------------------------------------

# these figs are simplified looks into the data and should not
# perfectly reflect the models we'll fit

for(foc in c("A", "B", "C")){
  foc_sp <- ifelse(foc == "A", "LB", 
                   ifelse(foc == "B", "AA", "CG"))
  temp <- fec_dat %>%
    filter(sp == foc & competitor %in% c("A", "ALONE")) %>%
    ggplot(aes(x = A_comp, y = fit)) + 
    geom_jitter(width = .25) + 
    stat_smooth(method = "nls",
                formula = y ~ a/(1+b*x),
                method.args = list(start = list(a = mean(fec_dat$fit[1:nrow(fec_dat)], na.rm = TRUE), b = .1)),
                se = FALSE, na.rm = TRUE, 
                color = colors[1]) + theme_classic() +
    xlab("Competitor Density") + 
    facet_grid(nitrogen ~ rhizobia, 
               labeller = labeller(nitrogen = c("nitrogen" = "Fertilized", 
                                                "zero nitrogen" = "Control"), 
                                   rhizobia = c("rhizobia" = "Inoculated", 
                                                "sterile" = "Uninoculated"))) + 
    ggtitle(paste(foc_sp, "LB", sep = "-")) + ylab("Fitness")
  
  ggsave(file.path("Figures", "Raw", "Fitness_Curves", 
                   paste(paste(foc_sp, "LB", sep = "-"), ".pdf", sep = "")), 
         plot = temp,
         width = 8, height = 6, units = "in")
  
  temp <- ggplot(data = fec_dat[which(fec_dat$sp == foc & fec_dat$competitor %in% c("B", "ALONE")),], 
                                    aes(x = B_comp, y = fit)) + 
    geom_jitter(width = .25) + 
    stat_smooth(method = "nls",
                formula = y ~ a/(1+b*x),
                method.args = list(start = list(a = mean(fec_dat$fit[1:nrow(fec_dat)], na.rm = TRUE), b = .1)),
                se = FALSE, na.rm = TRUE, 
                color = colors[1]) + theme_classic() +
    xlab("Competitor Density") + 
    facet_grid(nitrogen ~ rhizobia, 
               labeller = labeller(nitrogen = c("nitrogen" = "Fertilized", 
                                                "zero nitrogen" = "Control"), 
                                   rhizobia = c("rhizobia" = "Inoculated", 
                                                "sterile" = "Uninoculated"))) + 
    ggtitle(paste(foc_sp, "AA", sep = "-")) + ylab("Fitness")
  
  ggsave(file.path("Figures", "Raw", "Fitness_Curves", 
                   paste(paste(foc_sp, "AA", sep = "-"), ".pdf", sep = "")), 
         plot = temp,
         width = 8, height = 6, units = "in")
  
  
  temp <- ggplot(data = fec_dat[which(fec_dat$sp == foc & fec_dat$competitor %in% c("C", "ALONE")),], 
                                    aes(x = C_comp, y = fit)) + 
    geom_jitter(width = .25) + 
    stat_smooth(method = "nls",
                formula = y ~ a/(1+b*x),
                method.args = list(start = list(a = mean(fec_dat$fit[1:nrow(fec_dat)], na.rm = TRUE), b = .1)),
                se = FALSE, na.rm = TRUE, 
                color = colors[1]) + theme_classic() +
    xlab("Competitor Density") + 
    facet_grid(nitrogen ~ rhizobia, 
               labeller = labeller(nitrogen = c("nitrogen" = "Fertilized", 
                                                "zero nitrogen" = "Control"), 
                                   rhizobia = c("rhizobia" = "Inoculated", 
                                                "sterile" = "Uninoculated"))) + 
    ggtitle(paste(foc_sp, "CG", sep = "-")) + ylab("Fitness")
  
  ggsave(file.path("Figures", "Raw", "Fitness_Curves", 
                   paste(paste(foc_sp, "CG", sep = "-"), ".pdf", sep = "")), 
         plot = temp,
         width = 8, height = 6, units = "in")
}


# Mortality ---------------------------------------------------------------

fec_dat %>%
  group_by(nitrogen, rhizobia, sp2) %>%
  summarise(mort.rate = sum(fit == 0)/length(fit)) %>%
  ggplot(aes(x = nitrogen, y = mort.rate, color = rhizobia)) + 
  geom_point(position = position_dodge(width = 1), size = 2) + 
  facet_grid(~ sp2) + 
  scale_color_manual(values = colors, name = "Rhizobia", 
                     labels = c("Inoculated", "Uninoculated")) + 
  scale_x_discrete(labels = c("Fertilized", "Control"), name = "Nitrogen") + 
  theme_classic() + ylab("Mortality Proportion")

ggsave(file.path("Figures", "Raw", "Mortality_All.pdf"), 
       width = 8, height = 6, units = "in")

fec_dat %>%
  filter(density == "alone") %>%
  group_by(nitrogen, rhizobia, sp2) %>%
  summarise(mort.rate = sum(fit == 0)/length(fit)) %>%
  ggplot(aes(x = nitrogen, y = mort.rate, color = rhizobia)) + 
  geom_point(position = position_dodge(width = 1), size = 2) + 
  facet_grid(~ sp2) + 
  scale_color_manual(values = colors, name = "Rhizobia", 
                     labels = c("Inoculated", "Uninoculated")) + 
  scale_x_discrete(labels = c("Fertilized", "Control"), name = "Nitrogen") + 
  theme_classic() + ylab("Mortality Proportion")

ggsave(file.path("Figures", "Raw", "Mortality_Alone.pdf"), 
       width = 8, height = 6, units = "in")

