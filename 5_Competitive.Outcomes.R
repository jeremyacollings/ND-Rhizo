
########## Coexistence Results ###########

library(tidyverse)
library(tidybayes)
library(ggtern)

source("functions.R")

# Bring in data -----------------------------------------------------------

list2env(readRDS(file.path("RDS_Files", "derived_quantities.rds")), globalenv())

# Calculate Probabilities -------------------------------------------------

CX_df <- array2df(CX, "outcome")
winn_df <- array2df(winners, "winner")

full_df <- merge(CX_df, winn_df)

probs <- full_df %>%
  mutate(outcome = fct_recode(as.factor(outcome), 
                              coexistence = "1", 
                              exclusion = "2", 
                              priority_effects = "3")) %>%
  group_by(pair, rhizo, nitro) %>%
  summarise(pCX = sum(outcome == "coexistence")/n(), 
            pEX = sum(outcome == "exclusion")/n(), 
            pPE = sum(outcome == "priority_effects")/n())

probs

# Probability that rhizobia changes coexistence outcome? 

rhizo_switch <- full_df %>%
  mutate(outcome = fct_recode(as.factor(outcome), 
                              coexistence = "1", 
                              exclusion = "2", 
                              priority_effects = "3")) %>%
  replace_na(list(winner = 0)) %>%
  pivot_wider(id_cols = c(draw, pair, nitro), 
              names_from = c(rhizo), 
              values_from = c(outcome, winner)) %>%
  group_by(pair) %>%
  mutate(same_outcome = outcome_1 == outcome_2, 
         same_winner = winner_1 == winner_2) %>%
  summarise(p_switch_rhizo1 = sum(!same_outcome)/n(),
            p_switch_rhizo2 = sum(!same_winner, na.rm = TRUE)/n())

# Probability that nitrogen changes coexistence outcome? 

nitro_switch <- full_df %>%
  mutate(outcome = fct_recode(as.factor(outcome), 
                              coexistence = "1", 
                              exclusion = "2", 
                              priority_effects = "3")) %>%
  replace_na(list(winner = 0)) %>%
  pivot_wider(id_cols = c(draw, pair, rhizo), 
              names_from = c(nitro), 
              values_from = c(outcome, winner)) %>%
  group_by(pair) %>%
  mutate(same_outcome = outcome_1 == outcome_2, 
         same_winner = winner_1 == winner_2) %>%
  summarise(p_switch_nitro1 = sum(!same_outcome)/n(),
            p_switch_nitro2 = sum(!same_winner, na.rm = TRUE)/n())

# Probability that nitrogen changes whether or not rhizobia changes coexistence outcome? 

inter_switch <- full_df %>%
  mutate(outcome = fct_recode(as.factor(outcome), 
                              coexistence = "1", 
                              exclusion = "2", 
                              priority_effects = "3")) %>%
  replace_na(list(winner = 0)) %>%
  pivot_wider(id_cols = c(draw, pair), 
              names_from = c(rhizo, nitro), 
              values_from = c(outcome, winner)) %>%
  group_by(pair) %>%
  mutate(rhizo_same_outcome2 = outcome_2_2 == outcome_1_2, 
         rhizo_same_winner2 = winner_2_2 == winner_1_2, 
         rhizo_same_outcome1 = outcome_2_1 == outcome_1_1, 
         rhizo_same_winner1 = winner_2_1 == winner_1_1) %>%
  summarise(p_switch_inter1 = sum(rhizo_same_outcome2 != rhizo_same_outcome1)/n(),
            p_switch_inter2 = sum(rhizo_same_winner2 != rhizo_same_winner1)/n())

rhizo_switch %>%
  left_join(nitro_switch) %>%
  left_join(inter_switch)

# Graphs ------------------------------------------------------------------

full_dat <- array2df(ND, "ND") %>%
  left_join(array2df(FI, "FI")) %>%
  left_join(array2df(CX, "outcome")) %>%
  left_join(array2df(winners, "winner")) %>%
  mutate(treat = paste(rhizo, nitro))

summarized_dat <- full_dat %>%
  group_by(pair, treat) %>%
  median_qi(ND, FI, .width = 0.95) %>%
  # truncating values for graphing... keeping error bars within range
  mutate(FI.upper2 = ifelse(FI.upper > 75, 75, FI.upper), 
         ND.lower2 = ifelse(ND.lower < -5, -5, ND.lower))
  
samples <- sample(unique(full_dat$draw), 500)

# Good 'Ol Niche & Fitness Difference Coexistence Plot

full_dat %>%
  filter(draw %in% samples) %>%
  mutate(treat = paste(rhizo, nitro)) %>%
  ggplot(aes(x = ND, y = FI, shape = treat, color = pair)) + 
  geom_point(alpha = .1) + 
  geom_errorbar(data = summarized_dat, 
               aes(ymin = FI.lower, ymax = FI.upper2)) + 
  geom_errorbarh(data = summarized_dat, 
                aes(xmin = ND.lower2, xmax = ND.upper)) +
  stat_function(fun = function(x) -x+1, color = "black") + 
  stat_function(fun = function(x) 1/(-x+1), color = "black") + 
  geom_point(data = summarized_dat, size = 5) + 
  xlim(-6, 1) + ylim(1, 76) + 
  theme_classic(base_size = 15) + 
  scale_shape_manual(name = "Treatment", 
                     labels = c("Uninoculated, Control", 
                                "Uninoculated, Fertilized", 
                                "Inoculated, Control", 
                                "Inoculated, Fertilized"), 
                     values = c(1, 2, 16, 17)) + 
  scale_color_manual(name = "Pair", 
                     labels = c("AA.CG", "AA.LB", "CG.LB"), 
                     values = colors)

ggsave(file.path("Figures", "Coexistence_Metrics", "CX_graph.pdf"), 
       width = 8, height = 6, units = "in")

# Ternary Plot of Probabilities for Pair & Treatments

probs %>%
  mutate(treat = paste(rhizo, nitro)) %>%
  ggtern(aes(x = pCX, y = pEX, z = pPE,shape = treat, color = pair))  + 
  geom_point(size = 2) + 
  scale_shape_manual(name = "Treatment", 
                     labels = c("Uninoculated, Control", 
                                "Uninoculated, Fertilized", 
                                "Inoculated, Control", 
                                "Inoculated, Fertilized"), 
                     values = c(1, 2, 16, 17)) + 
  scale_color_manual(name = "Pair", 
                     labels = c("AA.CG", "AA.LB", "CG.LB"), 
                     values = colors) + 
  theme_classic() + 
  theme_zoom_T(.5) + 
  Tlab("EX") + Llab("CX") + Rlab("PE")

ggsave(file.path("Figures", "Coexistence_Metrics", "ternary.pdf"), 
       width = 8, height = 6, units = "in")
