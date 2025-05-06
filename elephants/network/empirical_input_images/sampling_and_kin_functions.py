import numpy as np
np.bool = np.bool_
from matplotlib import pyplot as plt
from scipy import stats
import sys
import itertools
import math
import pandas as pd
from PIL import Image, ImageDraw
rng = np.random.default_rng()

def find_pops(sample_array):
    '''
    Find and return all parent-offpsring pairs in sample_array
    '''
    # Find individuals in the sample that have a child in the sample
    ind_in_p1 = np.isin(sample_array.loc[:,'individual'].values, sample_array.loc[:,'parent1'])
    ind_in_p2 = np.isin(sample_array.loc[:,'individual'].values, sample_array.loc[:,'parent2'])

    # Mothers and Fathers
    mothers = sample_array[ind_in_p1]
    fathers = sample_array[ind_in_p2]
    parent_rows = pd.concat([mothers, fathers], axis = 0)
    parent_pairs = pd.DataFrame(columns = ['parent', 'year1', 'day1', 'x1', 'y1', 'child', 'year2', 'day2', 'x2', 'y2'])
    for i, parent in parent_rows.iterrows():
        parent_id = parent.loc['individual']
        children = sample_array[np.logical_or(sample_array.loc[:,'parent1'] == parent_id, sample_array.loc[:,'parent2'] == parent_id)]
        child_df = pd.DataFrame({'parent': [parent.loc['individual']]*len(children), 'year1': [parent.loc['year']]*len(children), 'day1': [parent.loc['day']]*len(children), 'x1': [parent.loc['x']]*len(children), 'y1': [parent.loc['y']]*len(children),
                                 'child': children.loc[:,'individual'], 'year2': children.loc[:,'year'], 'day2': children.loc[:,'day'], 'x2': children.loc[:,'x'], 'y2': children.loc[:,'y']})
        parent_pairs = pd.concat([parent_pairs, child_df], axis = 0)
    return(parent_pairs)

def find_recaptures(sample_array):
    '''
    Find all pairs of recaptures in sample_array
    '''
    recaptures = pd.DataFrame(columns = ['individual', 'year1', 'day1', 'x1', 'y1','year2', 'day2', 'x2', 'y2'])
    unique_inds = np.unique(sample_array.loc[:,'individual'], return_counts = True)
    recap_inds = unique_inds[0][unique_inds[1] > 1]
    recap_rows = sample_array[np.isin(sample_array.loc[:,'individual'], recap_inds)]
    for ind in recap_inds:
        ind_rows = recap_rows[recap_rows.loc[:,'individual'] == ind]
        # Sort rows by time
        ind_rows = ind_rows.sort_values(by = ['year', 'day'])
        # Add all pairs of recaptures
        for i in range(len(ind_rows)-1):
            for j in range(i+1, len(ind_rows)):
                row = ind_rows.iloc[i,:]
                row2 = ind_rows.iloc[j,:]
                new_recap_row = pd.DataFrame({'individual': [ind], 'year1': [row.loc['year']], 'day1': [row.loc['day']], 'x1': [row.loc['x']], 'y1': [row.loc['y']],
                                            'year2': [row2.loc['year']], 'day2': [row2.loc['day']], 'x2': [row2.loc['x']], 'y2': [row2.loc['y']]})
                recaptures = pd.concat([recaptures, new_recap_row], axis = 0)
    return(recaptures)

def find_sibs(sample_array):
    '''
    Find all pairs of half-siblings in sample_array
    '''
    # Find the parents of children in the sample and how many children in the sample they have (the parents don't need to be in the sample)
    # Mothers are parent1 and fathers are parent2
    all_mothers, all_mother_counts = np.unique(sample_array.loc[:,'parent1'].values, return_counts = True)
    all_fathers, all_father_counts = np.unique(sample_array.loc[:,'parent2'].values, return_counts = True)

    # Find parents with multiple children in the sample (parents of siblings)
    maternal_sib_parents = all_mothers[all_mother_counts > 1]
    paternal_sib_parents = all_fathers[all_father_counts > 1]
    sib_parents = np.concatenate([maternal_sib_parents, paternal_sib_parents])
    hs_pairs = pd.DataFrame(columns = ['hs1', 'year1', 'day1', 'x1', 'y1','hs2', 'year2', 'day2', 'x2', 'y2'])
    for parent_id in sib_parents:
        children = sample_array[np.logical_or(sample_array.loc[:,'parent1'] == parent_id,
                                            sample_array.loc[:,'parent2'] == parent_id)]
        children = children.sort_values(by = ['year', 'day'])
        # When one or both of the pair of half-sibs have been recaptured,
        # treat each recapture as a separate individual
        # and don't record pairs of the same individual
        for i in range(len(children)-1):
                for j in range(i+1, len(children)):
                    row1 = children.iloc[i,:]
                    row2 = children.iloc[j,:]
                    if(row1.loc['individual'] != row2.loc['individual']):
                        new_sib_row = pd.DataFrame({'hs1': [row1.loc['individual']], 'year1': [row1.loc['year']], 'day1': [row1.loc['day']], 'x1': [row1.loc['x']], 'y1': [row1.loc['y']],
                                                'hs2': [row2.loc['individual']],'year2': [row2.loc['year']], 'day2': [row2.loc['day']], 'x2': [row2.loc['x']], 'y2': [row2.loc['y']]})
                        hs_pairs = pd.concat([hs_pairs, new_sib_row], axis = 0)
    return(hs_pairs)

def draw_pairs(pairs, image, max_width, max_height, w, h, color = "white"):
    # pairs: An array of pairs to draw, with x1, y1, x2, y2 columns
    # image: An image object to add lines to
    # sample_array: Array containing information about each individual
    # max_width, max_height: Size of map used in SLiM simulation
    # w, h: Size of image
    for i, pair in pairs.iterrows():
        x1, y1 = pair[['x1', 'y1']].values
        x2, y2 = pair[['x2', 'y2']].values
        # Flip ys so origin is at lower left
        y1 = max_height - y1
        y2 = max_height - y2
        image.line([(x1*w/max_width, y1*h/max_height), (x2*w/max_width, y2*h/max_height)], fill =color, width = 0)
    