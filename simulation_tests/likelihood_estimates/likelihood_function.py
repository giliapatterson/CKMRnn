import numpy as np
from scipy import stats

# Survival, fecundity, and population size
def S(a, eta1_tilde, eta2_tilde, eta3_tilde, c):
    """
    Survival function
    """
    S = np.exp(-c*((eta1_tilde*a)**eta2_tilde + (eta1_tilde*a)**(1/eta2_tilde) + eta3_tilde - (eta1_tilde*(a-1))**eta2_tilde - (eta1_tilde*(a-1))**(1/eta2_tilde)))
    #S = np.exp(-c*((eta1_tilde*a)**eta2_tilde + (eta1_tilde*a)**(1/eta2_tilde) + eta3_tilde*a))
    return(S)
def f(a, nug1_tilde, nug2_tilde):
    """
    Fecundity function
    """
    return((1 + np.exp(-nug1_tilde*(a - nug2_tilde)))**(-1))
def N(params, max_a, max_t, c):
    """
    Number of males and females of ages 1, ..., max_a at years 1, ..., max_t
    """
    R0, eta1_tilde, eta2_tilde, eta3_tilde, nu11_tilde, nu12_tilde, nu21_tilde, nu22_tilde = params
    N_males = np.zeros((max_a + 1, max_t + 1))
    N_females = np.zeros((max_a + 1, max_t + 1))
    ages = np.arange(max_a + 1)
    years = np.arange(max_t + 1)
    # Set all time 0 population sizes to NA because time starts at 1
    N_males[:, 0] = None
    N_females[:, 0] = None
    # Set all age 0 population sizes to NA because age starts at 1
    N_males[0,:] = None
    N_females[0,:] = None
    # Initialize the population at time 1
    N_males[1, 1] = np.exp(R0)
    N_females[1, 1] = np.exp(R0)
    # For ages 2 to max_a, initialize as a stable age distribution
    for a in ages[2:]:
        S_amin1 = S(a - 1, eta1_tilde, eta2_tilde, eta3_tilde, c)
        N_males[a, 1] = N_males[a - 1, 1] * S_amin1
        N_females[a, 1] = N_females[a - 1, 1] * S_amin1
        #print("age: ", a, "time: 1", "Ns: ", N_males[a, 1], "   ", N_females[a, 1])
    for t in years[2:]:
        # Survival model
        for a in ages[2:]:
            S_amin1 = S(a - 1, eta1_tilde, eta2_tilde, eta3_tilde, c)
            N_males[a, t] = N_males[a - 1, t - 1] * S_amin1
            N_females[a, t] = N_females[a - 1, t - 1] * S_amin1
            #print("age: ", a, "time: ", t, "Ns: ", N_males[a, 1], "   ", N_females[a, 1])
        #print("age: ", 1, "time: ", t, "Ns: ", N_males[a, 1], "   ", N_females[a, 1])
        # Fecundity model
        n_new_offspring = np.sum([f(a, nu11_tilde, nu12_tilde) * N_females[a, t] for a in ages[2:,]])
        #print(n_new_offspring)
        N_males[1, t] = 0.5 * n_new_offspring
        N_females[1, t] = 0.5 * n_new_offspring
    return N_females, N_males

# Calculating probability of parent-offspring pair
def POP_probability(b1, y1, b2, g, N_females, N_males, nu11_tilde, nu12_tilde, nu21_tilde, nu22_tilde, max_a = 37):
    ages = np.arange(1,max_a + 1)
    if y1 >= b2 and b1 < y1 <= b1 + max_a and b1 <= b2 <= b1 + max_a:
        if g == 'F':
            nug1_tilde = nu11_tilde
            nug2_tilde = nu12_tilde
            N_g = N_females
        if g =='M':
            nug1_tilde = nu21_tilde
            nug2_tilde = nu22_tilde
            N_g = N_males
        # Total reproductive output in the offpsring birth year of all individuals of gender g alive at that time
        erro = np.sum([f(a, nug1_tilde, nug2_tilde) * N_g[a, b2] for a in ages[2:]])
        # Fecundity of parent in offspring birth year
        f_parent = f(b2 - b1 + 1, nug1_tilde, nug2_tilde)
        p_pop = f_parent/erro
    else:
        p_pop = 0
    return(p_pop)

# Calculating probability of half-sibling pair
def HS_probability(b1, b2, g, N_females, N_males, nu11_tilde, nu12_tilde, nu21_tilde, nu22_tilde, eta1_tilde, eta2_tilde, eta3_tilde, c, max_a = 37):
    if b2 < b1:
        return(0)
    ages = np.arange(1,max_a + 1)
    if g == 'F':
        nug1_tilde = nu11_tilde
        nug2_tilde = nu12_tilde
        N_g = N_females
    if g =='M':
        nug1_tilde = nu21_tilde
        nug2_tilde = nu22_tilde
        N_g = N_males
    p_hs = 0
    # Sum over all possible ages
    for a in ages:
        # Total number of individuals of sex g and age a in the older sibling birth year
        N_ab1g = N_g[a, b1]
        # Total expected reproductive output of all individuals of sex g in the older sibling birth year
        erro_b1 = np.sum([f(a_prime, nug1_tilde, nug2_tilde) * N_g[a_prime, b1] for a_prime in ages[2:]])
        # Probability that an individual of age a in the older sibling birth year is the parent of the older sibling
        p_parent_older_sib = f(a, nug1_tilde, nug2_tilde)/erro_b1
        # Probability that an individual of age a in the older sibling birth year survived until the younger sibling birth year
        p_survived = np.prod([S(a + y - b1, eta1_tilde, eta2_tilde, eta3_tilde, c) for y in np.arange(b1, b2)])
        # Total expected reproductive output of all individuals of sex g in the younger sibling birth year
        erro_b2 = np.sum([f(a_prime, nug1_tilde, nug2_tilde) * N_g[a_prime, b2] for a_prime in ages[2:]])
        # Probability that an individual of age a in the older sibling birth year is the parent of the younber sibling
        p_parent_younger_sib = f(a + b2 - b1, nug1_tilde, nug2_tilde)/erro_b2
        p_hs += N_ab1g * p_parent_older_sib * p_survived * p_parent_younger_sib
    return(p_hs)

def nll(params, survival_prior_params, fecundity_prior_params, parent_counts, sibling_counts, parent_comparisons, sibling_comparisons, c = 1.111, max_a = 37, max_t = 60):
    """
    Log likelihood function of params given parent and sibling pair counts
    Parameters to estimate: R0, eta1_tilde, eta2_tilde, eta3_tilde, nu11_tilde, nu12_tilde, nu21_tilde, nu22_tilde
    Observations: parent_counts, sibling_counts, parent_comparisons, sibling_comparisons
    """
    # Unpack the parameters
    R0, eta1_tilde, eta2_tilde, eta3_tilde, nu11_tilde, nu12_tilde, nu21_tilde, nu22_tilde = params
    eta1, eta2, eta3, sigma1, sigma2, sigma3 = survival_prior_params
    nu11, nu12, nu21, nu22, tau11, tau12, tau21, tau22 = fecundity_prior_params
    # Population size model
    N_females, N_males = N(params, max_a, max_t, c)
    # Parent section of the negative log likelihood
    PO_nll = 0
    for b1, y1, b2, g in parent_comparisons.keys():
        p_pop = POP_probability(b1, y1, b2, g, N_females = N_females, N_males = N_males, nu11_tilde = nu11_tilde, nu12_tilde = nu12_tilde, nu21_tilde = nu21_tilde, nu22_tilde = nu22_tilde, max_a = max_a)
        n = parent_comparisons[(b1, y1, b2, g)]
        m = parent_counts.get((b1, y1, b2, g), 0)
        if p_pop != 0:
            PO_nll += -m*np.log(n*p_pop) + n*p_pop
    #print(PO_nll)
    # Sibling section of the negative log likelihood
    HS_nll = 0
    for b1, b2, g in sibling_comparisons.keys():
        p_hs = HS_probability(b1, b2, g, N_females = N_females, N_males = N_males, nu11_tilde = nu11_tilde, nu12_tilde = nu12_tilde, nu21_tilde = nu21_tilde, nu22_tilde = nu22_tilde, eta1_tilde = eta1_tilde, eta2_tilde = eta2_tilde, eta3_tilde = eta3_tilde, c = c, max_a = max_a)
        n_prime = sibling_comparisons[(b1, b2, g)]
        m_prime = sibling_counts.get((b1, b2, g), 0)
        if p_hs != 0:
            HS_nll += -m_prime*np.log(n_prime*p_hs) + n_prime*p_hs
    #print(HS_nll)
    # Negative log likelihood for priors on the survival parameters
    survival_nll = -(stats.norm.logpdf(eta1_tilde, eta1, sigma1) + stats.norm.logpdf(eta2_tilde, eta2, sigma2) + stats.norm.logpdf(eta3_tilde, eta3, sigma3))
    # Negative log likelihood for priors on the fecundity parameters
    fecundity_nll = -(stats.norm.logpdf(nu11_tilde, nu11, tau11) + stats.norm.logpdf(nu12_tilde, nu12, tau12) +
                    stats.norm.logpdf(nu21_tilde, nu21, tau21) + stats.norm.logpdf(nu22_tilde, nu22, tau22))
    nll = PO_nll + HS_nll + survival_nll + fecundity_nll
    return(nll)