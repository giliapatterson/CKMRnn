library(tidyverse)

# Load test simulations
test_sims <- read_csv("ci_calibration_results/test_sims_ranks.csv")

# Make a qq plot of ranks compared to uniform distribution
ggplot(test_sims, aes(sample = rank_normalized)) +
  stat_qq(distribution = qunif) +
  ylab("Rank of true value in bootstrap replicates") +
  xlab("Theoretical quantiles of uniform distribution") +
  theme_bw(base_size = 14) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
ggsave("../paper_figures/ci_calibration_qqplot.pdf", width = 6, height = 6)

ggplot(test_sims, aes(x = rank_normalized)) +
  geom_histogram(stat = "density") +
  xlab("Rank of true value in bootstrap replicates") +
  ylab("Density") +
  theme_minimal()

mean(test_sims$in_95_CI)
mean(test_sims$in_50_CI)
