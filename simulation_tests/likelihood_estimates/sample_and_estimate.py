import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from scipy import optimize
import pandas as pd
import nlopt
import argparse

# Functions for sampling and finding kin
import sampling_and_kin_functions as skf
import likelihood_function as lf

def sample_and_estimate(parents_file, bias, noptsteps):
    parents = pd.read_csv(parents_file)
    # True values of parameters for SLiM simulation
    max_a = 37
    max_t = 60
    R0 = 8
    eta1_tilde = np.exp(-2.904)
    eta2_tilde = 1 + np.exp(0.586)
    eta3_tilde = np.exp(-2.579)
    c = 1.111
    nu11_tilde = 1.264
    nu12_tilde = 5.424
    nu21_tilde = 1.868
    nu22_tilde = 6.5
    true_params = [R0, eta1_tilde, eta2_tilde, eta3_tilde, nu11_tilde, nu12_tilde, nu21_tilde, nu22_tilde]

    # Add birth year to the data frame
    # Because individuals are 1 in the year they are born, add 1
    parents["birth_year"] = parents["sampling_time"] - parents["age"] + 1

    # Define sampling intensity grid with one side [bias] times as likely to be sampled as the other
    sampling_intensity =  np.repeat([np.linspace(1, bias, 10)], 10, axis = 0)

    # Sample according to the grid and return realized sampling size and sampled individuals
    ss, sample_rows = skf.sample_grid(parents, sampling_intensity, 2000, 10, 10)
    sample_parents = parents.iloc[sample_rows]

    # Find individuals in the sample that have a child in the sample
    ind_in_p1 = np.isin(sample_parents.loc[:,'individual'].values, sample_parents.loc[:,'parent1'])
    ind_in_p2 = np.isin(sample_parents.loc[:,'individual'].values, sample_parents.loc[:,'parent2'])

    # Mothers and Fathers
    mothers = sample_parents[ind_in_p1]
    fathers = sample_parents[ind_in_p2]
    mother_list = mothers.loc[:,'individual'].values
    father_list = fathers.loc[:,'individual'].values

    # Find the parents of children in the sample and how many children in the sample they have (the parents don't need to be in the sample)
    # Mothers are parent1 and fathers are parent2
    all_mothers, all_mother_counts = np.unique(sample_parents.loc[:,'parent1'].values, return_counts = True)
    all_fathers, all_father_counts = np.unique(sample_parents.loc[:,'parent2'].values, return_counts = True)

    # Find parents with multiple children in the sample (parents of siblings)
    maternal_sib_parents = all_mothers[all_mother_counts > 1]
    paternal_sib_parents = all_fathers[all_father_counts > 1]

    # Record POPs and half-sibling pairs
    # POPs are [parent, offspring], sibs are [sib1, sib2] where sib1 is the sibling with a lower id number (which might mean younger but I'd need to check)
    maternal_pops = skf.find_POPs_or_sibs(mother_list, sample_parents, "PO")
    paternal_pops = skf.find_POPs_or_sibs(father_list, sample_parents, "PO")
    maternal_sibs = skf.find_POPs_or_sibs(maternal_sib_parents, sample_parents, "HS")
    paternal_sibs = skf.find_POPs_or_sibs(paternal_sib_parents, sample_parents, "HS")

    # Full sibling pairs appear in both maternal and paternal arrays
    if np.any([x in paternal_sibs for x in maternal_sibs]):
        full_sibs = maternal_sibs[np.where([x in paternal_sibs for x in maternal_sibs])[0][0]]
    else:
        full_sibs = []

    # Initialize empty dictionaries to keep track of counts
    parent_counts = {}
    sibling_counts = {}

    # For each parent birth year b1, parent capture year y1, offspring birth year b2, and parent sex g
    # Count the number of parent offspring pairs
    pairs = np.concatenate([maternal_pops, paternal_pops])
    sample_array = sample_parents
    # pairs: An array of tuples of individual ids
    # sample_array: Array containing information about each individual
    for pair in pairs:
        index1 = np.where(sample_array.loc[:,'individual'] == pair[0])
        index2 = np.where(sample_array.loc[:,'individual'] == pair[1])
        row1 = sample_array.iloc[index1]
        row2 = sample_array.iloc[index2]
        b1 = row1['birth_year'].values[0]
        b2 = row2['birth_year'].values[0]
        y1 = row1['sampling_time'].values[0]
        g =  row1['sex'].values[0]
        parent_counts[(b1, y1, b2, g)] = parent_counts.get((b1, y1, b2, g), 0) + 1

    # For each older sibling birth year b1, younger sibling birth year b2, and parent sex g,
    # count the number of sibling pairs
    pairs = maternal_sibs
    # pairs: An array of tuples of individual ids
    for pair in pairs:
        index1 = np.where(sample_array.loc[:,'individual'] == pair[0])
        index2 = np.where(sample_array.loc[:,'individual'] == pair[1])
        row1 = sample_array.iloc[index1]
        row2 = sample_array.iloc[index2]
        sib1_b = row1['birth_year'].values[0]
        sib2_b = row2['birth_year'].values[0]
        b1, b2 = sorted([sib1_b, sib2_b])
        g = 'F'
        sibling_counts[(b1, b2, g)] = sibling_counts.get((b1, b2, g), 0) + 1
    pairs = paternal_sibs
    # pairs: An array of tuples of individual ids
    for pair in pairs:
        index1 = np.where(sample_array.loc[:,'individual'] == pair[0])
        index2 = np.where(sample_array.loc[:,'individual'] == pair[1])
        row1 = sample_array.iloc[index1]
        row2 = sample_array.iloc[index2]
        sib1_b = row1['birth_year'].values[0]
        sib2_b = row2['birth_year'].values[0]
        b1, b2 = sorted([sib1_b, sib2_b])
        g = 'M'
        sibling_counts[(b1, b2, g)] = sibling_counts.get((b1, b2, g), 0) + 1

    # Number of individuals born in each year
    b_grouped = sample_parents.groupby(['birth_year'])
    b_counts = {key: len(value) for key, value in b_grouped.groups.items()}
    # Number of individuals of each birth year, capture year, and sex
    byg_grouped = sample_parents.groupby(['birth_year', 'sampling_time', 'sex'])
    byg_counts = {key: len(value) for key, value in byg_grouped.groups.items()}

    # Number of parent offspring pair comparisons
    # Initialize dictionary
    parent_comparisons = {}
    max_a = 37
    # Loop through all parents
    for b1, y1, g in byg_counts.keys():
        # Loop through all offpsring
        for b2 in b_counts.keys():
            if (y1 >= b2 and b1 < y1 <= b1 + max_a and b1 <= b2 <= b1 + max_a):
                parent_comparisons[(b1, y1, b2, g)] = byg_counts[(b1, y1, g)] * b_counts[b2]                   
    for b1, y1, b2, g in parent_counts.keys():
        if(parent_comparisons[(b1, y1, b2, g)] < parent_counts[(b1, y1, b2, g)]):
            print("Error, more PO pairs than possible")
            break

    # Number of half-sibling pair comparisons
    # Initialize dictionary
    sibling_comparisons = {}
    # Loop through all older siblings
    for b1 in b_counts.keys():
        # Loop through all younger siblings
        for b2 in b_counts.keys():
            if b2 >= b1:
                sibling_comparisons[(b1, b2, 'F')] = b_counts[b1] * b_counts[b2] 
                sibling_comparisons[(b1, b2, 'M')] = b_counts[b1] * b_counts[b2]                  
    for b1, b2, g in sibling_counts.keys():
        if(sibling_comparisons[(b1, b2, g)] < sibling_counts[(b1, b2, g)]):
            print("Error, more sibling pairs than possible")
            break
    
    # Define the objective function
    def nll_function(x, grad):
        """
        x = R0
        grad is empty
        """
        survival_priors = 0.055, 2.8, 0.076, np.exp(0.15), np.exp(0.25), np.exp(0.5)
        fecundity_priors = 1.264, 5.424, 1.868, 6.5, 0.4*1.264, 0.4*5.424, 0.4*1.868, 0.4*6.5
        eta1, eta2, eta3, sigma1, sigma2, sigma3 = survival_priors
        nu11, nu12, nu21, nu22, tau11, tau12, tau21, tau22 = fecundity_priors
        nll_val = lf.nll([x, eta1, eta2, eta3, nu11, nu12, nu21, nu22], survival_priors, fecundity_priors, parent_counts, sibling_counts, parent_comparisons, sibling_comparisons)
        print(nll_val)
        return nll_val

    opt = nlopt.opt(nlopt.LN_BOBYQA, 1)
    opt.set_min_objective(nll_function)
    # Parameters: R0, eta1_tilde, eta2_tilde, eta3_tilde, nu11_tilde, nu12_tilde, nu21_tilde, nu22_tilde
    # R0 > 0, all other parameters are unconstrained
    opt.set_lower_bounds(7)
    opt.set_upper_bounds(12)
    opt.set_maxeval(noptsteps) # Set max iterations, print likelihood to see that it's decreasing
    #print(nll_function([10, eta1, eta2, eta3, nu11, nu12, nu21, nu22]))
    x = opt.optimize([10])
    print(x)

    estimatedR0 = x[0]
    estimatedN = np.exp(x[0])
    return(estimatedR0, estimatedN)