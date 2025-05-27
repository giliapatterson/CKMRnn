library(tidyverse)
theme_set(theme_gray(base_size = 22))
# Constants and survival probability from Conn et al. code
# https://github.com/pconn/CKMR/blob/master/R/simulate_spatial.R
a.haz = exp(-2.904)
b.haz = 1 + exp(0.586)
c.haz = exp(-2.579)
phi_RAW = function(a.haz,b.haz,c.haz,haz.mult,Age){
  exp(-haz.mult*((a.haz*Age)^b.haz + (a.haz*Age)^(1/b.haz) + c.haz - (a.haz*(Age-1))^b.haz - (a.haz*(Age-1))^(1/b.haz) ))
}

# Plot survival, breeding, and mating probabilities
ages = seq(1, 38, 0.1)
p_breed = 1/((1 + exp(-1.264*(ages-5.424))))
S_conn =  phi_RAW(a.haz, b.haz, c.haz, 1.111, ages)
p_mate = 1/((1+exp(-1.868*(ages - 6.5))))

parameters <- data.frame(age = ages, p_breed = p_breed, p_mate = p_mate, S_conn = S_conn)
parameters <- pivot_longer(parameters, cols = c(p_breed, p_mate, S_conn))
ggplot(parameters, aes(x = age, y = value, color = name)) +
  geom_line() +
  ylim(0, 1) +
  scale_color_discrete(name = "Parameter", labels = c("Breeding", "Mating", "Survival")) +
  xlab("Age") +
  ylab("Probability")

## Figure out the multiplier that gives constant population size(C in the paper, haz.mult in Conn code) ##

# Function to make transition matrix for Leslie matrix model for females
transition_matrix = function(a.haz, b.haz, c.haz, C, yearly_ages){
  p_breed = 1/((1 + exp(-1.264*(yearly_ages-5.424))))
  S_conn =  phi_RAW(a.haz, b.haz, c.haz, C, yearly_ages)
  p_mate = 1/((1+exp(-1.868*(yearly_ages - 6.5))))
  # Survival
  Surv = diag(S_conn)
  # First row is female breeding probability divided by 2
  Fec = matrix(p_breed/2, nrow = 1)
  # Leslie matrix
  M = rbind(Fec, Surv[1:(nrow(Surv)-1),])
  return(M)
}

# Leslie matrix for C = 1
yearly_ages = 1:38
a.haz = exp(-2.904)
b.haz = 1 + exp(0.586)
c.haz = exp(-2.579)
M = transition_matrix(a.haz, b.haz, c.haz, 1, yearly_ages)
# Find eigenvalues
eigen(M)$values[1]

# Which value of C gives a first eigenvalue of 1?

# Values to try
C_vals <- seq(1, 1.3, 0.001)
# First eigenvalue (want it to be one)
first_eigenvalue <- rep(NA, length(C_vals))
# Possible ages consistent with the Conn paper
yearly_ages = 1:38
# Parameters
a.haz = exp(-2.904)
b.haz = 1 + exp(0.586)
c.haz = exp(-2.579)
i = 1
for (C in C_vals){
  M_C = transition_matrix(a.haz, b.haz, c.haz, C, yearly_ages)
  first_eigenvalue[i] = eigen(M_C)$values[1]
  i = i + 1
}
plot(C_vals, first_eigenvalue)

# Which C gives first eigenvalue closest to 1?
C_i = which.min(abs(first_eigenvalue - 1))
C_1 = C_vals[C_i]

# Check that the eigenvalue for C_1 is close to 1
M_1 = transition_matrix(a.haz, b.haz, c.haz, C_1, yearly_ages)
print(eigen(M_1)$values[1]) # Should be close to 1

print(M_1)

# Get stable age distribution for C_1
stable_age_distribution = as.numeric(eigen(M_1)$vectors[,1]/sum(eigen(M_1)$vectors[,1]))
dput(round(stable_age_distribution, 3))
plot(yearly_ages, stable_age_distribution)
