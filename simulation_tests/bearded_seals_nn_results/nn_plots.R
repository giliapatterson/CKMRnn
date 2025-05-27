library(tidyverse)
library(gridExtra)
theme_set(theme_gray(base_size = 50))

# Results for less than or greater than this size
breakpoint = 12000
# Results for bearded seal network trained and tested on constant population size
con_con <- read.csv("constant_size_network.csv")
con_con <- con_con %>% mutate(breakpoint = ifelse(truth < breakpoint, 'N<12000', 'N>12000'))
# Results for bearded seal network trained on constant population size and tested on changing population size
con_change <- read.csv("misspecified_trend.csv")
con_change <- mutate(con_change, surv_mult = factor(surv_mult)) %>%
  mutate(breakpoint = ifelse(true_N_final < breakpoint, 'N<12000', 'N>12000'))
# Results for bearded seal network both trained and tested on changing population size
change_change <- read.csv("varying_trend_network.csv")
change_change <- mutate(change_change, surv_mult = factor(surv_mult))
change_change <- change_change %>%
  mutate(breakpoint = ifelse(true_N_final < breakpoint, 'N<12000', 'N>12000'))
# Compute MSE and mean relative error (as percent of truth)
con_con <- mutate(con_con, error = bearded_nn_pred - truth,
                  rel_error = error/truth, bias = factor(bias))
con_change<- mutate(con_change, error = bearded_nn_pred - true_N_final,
                    rel_error = error/true_N_final, bias = factor(bias))
change_change<- mutate(change_change, error = bearded_nn_pred - true_N_final,
                       rel_error = error/true_N_final, bias = factor(bias))

### Results for constant population size ###
# True vs predicted
bias_labels <- c('1' = 'Random', '16.5' = 'Medium bias', '31.875' = 'High bias')
tvsp <- ggplot(con_con, aes(x = truth, y = bearded_nn_pred)) +
  geom_point() +
  xlab("True N") +
  ylab("Predicted N") +
  geom_abline(slope = 1) +
  facet_wrap(~bias, labeller = as_labeller(bias_labels))
violin_df <- con_con
violin_df$breakpoint = 'All'
violin_df <- bind_rows(con_con, violin_df)
violin <- ggplot(violin_df, aes(x = breakpoint, y = rel_error)) +
  geom_violin() +
  geom_hline(yintercept = 0) +
  facet_wrap(~bias, labeller = as_labeller(bias_labels)) +
  xlab("Test set population sizes") +
  ylab("Relative error")
g <- arrangeGrob(tvsp, violin)
ggsave("constant_constant.png", g, width = 30, height = 20)
ggsave("constant_constant_tvsp.png", tvsp, width = 30, height = 10)
ggsave("constant_constant_error.png", violin, width = 30, height = 10)

# Relative error
ggplot(con_con, aes(x = truth, y = rel_error)) +
  geom_point() +
  xlab("True N") +
  ylab("Error") +
  geom_abline(slope = 1) +
  facet_wrap(~bias, labeller = "label_both") +
  geom_smooth()

## Results for constant-size training, changing size test
# Filter to use only bias of 15.5 and only increasing and decreasing trends
con_change15 <- filter(con_change, bias == 15.5, surv_mult != 1)
# Labels
surv_labels <- c('0.99' = "Decreasing", '1' = "Constant", "1.01" = "Increasing")
# True vs predicted
tvsp <- ggplot(con_change15, aes(x = true_N_final, y = bearded_nn_pred)) +
  geom_point() +
  xlab("True Final N") +
  ylab("Predicted N") +
  geom_abline(slope = 1) +
  facet_wrap(~surv_mult, scales = "free_y", labeller = as_labeller(surv_labels))
# Relative error
violin <- ggplot(con_change, aes(x = surv_mult, y = rel_error)) +
  geom_violin() +
  scale_x_discrete(labels = as_labeller(surv_labels)) +
  geom_hline(yintercept = 0) +
  xlab("Population trend") +
  ylab("Relative error")
g <- arrangeGrob(tvsp, violin, nrow = 1, widths = c(2/3, 1/3))
ggsave("constant_changing.png", g, width = 30, height = 10)
ggsave("constant_changing_tvsp.png", tvsp, width = 20, height = 10)
ggsave("constant_changing_error.png", violin, width = 20, height = 10)

## Results for changing-size training, changing size test
# Filter to use only bias of 15.5
change_change15 <- filter(change_change, bias == 15.5)
# Labels
surv_labels <- c('0.99' = "Decreasing", '1' = "Constant", "1.01" = "Increasing")
# True vs predicted
tvsp <- ggplot(change_change15, aes(x = true_N_final, y = bearded_nn_pred)) +
  geom_point() +
  xlab("True Final N") +
  ylab("Predicted N") +
  geom_abline(slope = 1) +
  facet_wrap(~surv_mult, scales = "free_y", labeller = as_labeller(surv_labels))
# Relative error
violin <- ggplot(change_change, aes(x = surv_mult, y = rel_error)) +
  geom_violin() +
  scale_x_discrete(labels = as_labeller(surv_labels)) +
  geom_hline(yintercept = 0) +
  xlab("Population trend") +
  ylab("Relative error")
g <- arrangeGrob(tvsp, violin, nrow = 1, widths = c(3/4, 1/4))
ggsave("changing_changing.png", g, width = 40, height = 10)
ggsave("changing_changing_tvsp.png", tvsp, width = 30, height = 10)
ggsave("changing_changing_error.png", violin, width = 20, height = 10)

## Table of mean absolute relative error ##
# Columns: Training set, test set bias, test set trend, mrae overall, mrae < 12000, mrae > 12000
# Constant-constant
con_mrae <- con_con %>% group_by(bias) %>% summarize(mrae = mean(abs(rel_error)),
                                                     mre = mean(rel_error))
con_mrae <- con_con %>% group_by(bias, breakpoint) %>%
  summarize(mrae = mean(abs(rel_error)), mre = mean(rel_error)) %>%
  pivot_wider(names_from = breakpoint, values_from = c(mrae, mre)) %>%
  full_join(con_mrae)  %>%
  mutate(training = "constant", surv_mult = 1)
# Constant-changing
con_change_mrae <- con_change %>% filter(surv_mult !=1) %>% group_by(bias, surv_mult) %>%
  summarize(mrae = mean(abs(rel_error)), mre = mean(rel_error)) %>%
  mutate(training = 'constant')
# Changing-changing
change_change_mrae <- change_change %>% group_by(bias, surv_mult) %>%
  summarize(mrae = mean(abs(rel_error)), mre = mean(rel_error)) %>%
  mutate(training = 'varying')

# Table for constant-constant
con_con_table <- con_mrae %>%
  ungroup() %>%
  mutate(trend = case_when(surv_mult == 0.99 ~ 'decreasing',
                           surv_mult == 1.01 ~ 'increasing',
                           surv_mult == 1 ~ 'constant',
                           .default = "constant"),
         bias = as.numeric(levels(bias))[bias],
         spatial_bias = case_when(bias ==1 ~ 'random',
                          bias > 10 & bias < 20 ~ 'medium',
                          bias > 20 ~ 'high')) %>%
  select(spatial_bias, mrae, 'mrae_N<12000', 'mrae_N>12000',
         mre, 'mre_N<12000', 'mre_N>12000')
knitr::kable(con_con_table, format = 'latex', digits = 3,
             col.names = c("Spatial sampling bias",
                           "Mean relative absolute error",
                           "MRAE N<12000",
                           "MRAE N>=2000",
                           "Mean relative error",
                           "MRE N<12000",
                           "MRE N>=2000"))

# Table for constant-changing
con_change_table <- con_change_mrae %>%
  ungroup() %>%
  mutate(trend = case_when(surv_mult == 0.99 ~ 'decreasing',
                           surv_mult == 1.01 ~ 'increasing',
                           surv_mult == 1 ~ 'constant',
                           .default = "constant"),
         bias = as.numeric(levels(bias))[bias],
         spatial_bias = case_when(bias ==1 ~ 'random',
                                  bias > 10 & bias < 20 ~ 'medium',
                                  bias > 20 ~ 'high')) %>%
  select(trend, spatial_bias, mrae, mre)
knitr::kable(con_change_table, format = 'latex', digits = 3,
             col.names = c("Population trend",
                           "Spatial sampling bias",
                           "Mean relative absolute error",
                           "Mean relative error"))

# Table for changing-changing
change_change_table <- change_change_mrae %>%
  ungroup() %>%
  mutate(trend = case_when(surv_mult == 0.99 ~ 'decreasing',
                           surv_mult == 1.01 ~ 'increasing',
                           surv_mult == 1 ~ 'constant',
                           .default = "constant"),
         bias = as.numeric(levels(bias))[bias],
         spatial_bias = case_when(bias ==1 ~ 'random',
                                  bias > 10 & bias < 20 ~ 'medium',
                                  bias > 20 ~ 'high')) %>%
  select(trend, spatial_bias, mrae, mre)
knitr::kable(change_change_table, format = 'latex', digits = 3,
             col.names = c("Population trend",
                           "Spatial sampling bias",
                           "Mean relative absolute error",
                           "Mean relative error"))


