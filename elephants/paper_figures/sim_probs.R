library(tidyverse)

# Parameters
# Fixed parameters
FAGEREPRO = 20
# Breeding probabilities of females
age = 0:60
last_rep = -6:0
tick = 0
gap = tick - last_rep
females = expand_grid(age, gap) %>%
  mutate(page =  1/(1+exp(-0.5*(age-FAGEREPRO))),
         pgap = 1/(1+exp(-2*(gap-5))),
         pbreed = pmin(1, page*pgap))

ggplot(females, aes(x = age, y = page)) +
  geom_point() +
  geom_line()
ggplot(females, aes(x = gap, y = pgap)) +
  geom_point() +
  geom_line()
ggplot(females, aes(x = age, y = pbreed, color = as.factor(gap))) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d() +
  theme_grey(base_size = 30) +
  xlab("Age (years)") +
  ylab("Probability of fertility") +
  labs(color = "Years since\n last reproduction")
ggsave("female_repro_probs.png", width = 10, height = 10)

females %>% group_by(gap) %>% summarize(max = max(pbreed))

### Survival probabilities
survival_probs <- read_csv("../one_year_survival.csv") %>%
  pivot_longer(cols = c('F', 'M'), names_to = "Sex", values_to = "Probability")
ggplot(survival_probs, aes(x = Age, y = Probability, color = Sex)) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d() +
  theme_grey(base_size = 30) +
  xlab("Age (years)") +
  ylab("Probability of survival") +
  labs(color = "Sex")

table <-survival_probs %>% mutate(age_class = floor(survival_probs$Age/5)) %>%
  group_by(age_class, Sex) %>% summarize(prob = unique(Probability),
                                         name = paste(min(Age), max(Age), sep="-")) %>%
  pivot_wider(names_from = Sex, values_from = prob) %>% ungroup %>%
  select(name, 'F', 'M') %>%
  filter(name != '60-60') %>%
  rename(Age = name) %>%
  mutate(Age = if_else(Age == "55-59", "\\geq 55", Age))

knitr::kable(table, col.names = c("Age", "Female", "Male"), digits = 3,
             format = 'latex', escape = FALSE)

