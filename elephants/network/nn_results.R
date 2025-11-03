library(tidyverse)
library(cowplot)
theme_set(theme_gray(base_size = 14))

## Training history
hist <- read.csv("network_output/training_hist.csv")
hist <- hist %>% rename(Epoch = train_x, Loss = train_loss)

ggplot(hist, aes(x = Epoch, y = Loss)) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d() +
  xlab("Training Epoch") +
  ylab("Loss (Mean Squared Error)")
ggsave("network_output/training_hist.png")

## Performance on held out test set of simulations ##
test_set <- read.csv("network_output/test_res.csv")
test_set <- mutate(test_set, error = pred - N,
                  rel_error = error/N)

tvsp <- ggplot(test_set, aes(x = N, y = pred)) +
  geom_point(size = 0.5, stroke = 0.2) +
  xlab("True N") +
  ylab("Predicted N") +
  geom_abline(slope = 1, linewidth = 0.2)

violin <- ggplot(test_set, aes(x = "", y = rel_error)) +
  geom_violin(linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  ylab("Relative error") +
  xlab("Test set")

g <- plot_grid(tvsp, violin, align = "hv", ncol = 2,
               labels = c("(a)", "(b)"),
               label_x = 0, label_y = 0,
               hjust = -0.5, vjust = -0.5)
ggsave("../paper_figures/figure6.pdf", g, width = 6, height = 3)


# Mean relative absolute error
(elephant_mrae <- mean(abs(test_set$rel_error)))

## Parametric bootstrapping
bootstrap <- read.csv("network_output/bootstrap_reps.csv")
# Point estimate
estimate <- unique(bootstrap$estimate)[1]
# 95% Confidence interval
ci = quantile(bootstrap$pred, probs = c(0.025, 0.975))
lower = ci[1]
upper = ci[2]
labels = c(paste("95% CI Lower Bound:", round(lower, 0)),
           paste("95% CI Upper Bound:", round(upper, 0)),
           paste("Estimate:", round(estimate, 0)))

ggplot(bootstrap, aes(x = pred)) +
  geom_histogram(binwidth = 20) +
  geom_vline(xintercept = estimate, linewidth = 0.2) +
  geom_vline(xintercept = lower, linewidth = 0.2) +
  geom_vline(xintercept = upper, linewidth = 0.2) +
  xlab("Predicted N") +
  ylab("Count") +
  annotate(geom = "text",
          label = labels,
          x = c(lower, upper, estimate),
          y = 30,
          angle = 90,
          vjust = -0.5,
          size = 3)
ggsave("../paper_figures/figure7.pdf", width = 4, height = 4)




