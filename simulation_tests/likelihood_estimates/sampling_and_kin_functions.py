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

# Sample from the individuals according to a sampling intensity grid. 
# Individuals from all years are lumped together but I can change this later

def sample_cell(sample_parents, xmin, xmax, ymin, ymax, nmax):
    # Sample from individuals within one grid cell
    in_i = np.where(np.logical_and(np.logical_and(sample_parents.loc[:,'x'] <= xmax, sample_parents.loc[:,'x'] >= xmin),
                np.logical_and(sample_parents.loc[:,'y'] >= ymin, sample_parents.loc[:,'y'] <= ymax)))[0]
    # Sample size is nmax or the number of individuals in the area
    ss = min(len(in_i), nmax)
    sample_rows = rng.choice(in_i, ss, replace = False)
    return(ss, sample_rows)

def sample_grid(individuals, sampling_intensity, n, width = 10, height = 10):
    # Sample from each grid cell
    # Sample size 
    ss = np.array(np.ceil(n*sampling_intensity/np.sum(sampling_intensity)), dtype = int)
    # Keep track of actual sample sizes (since if there are fewer than ss individuals in a grid cell, all are sampled)
    realized_ss = np.zeros(ss.shape, dtype = int)
    
    # Width and height of grid cells
    x_cells = ss.shape[0]
    y_cells = ss.shape[1]
    dx = width/x_cells
    dy = height/y_cells
    
    # Keep track of sampled rows
    #sample_rows = np.empty(round(np.sum(ss)), dtype = int)
    sample_rows = np.empty(0, dtype = int)
    # Sample from each grid cell
    for ix, iy in np.ndindex(ss.shape):
        xmin = dx*ix
        xmax = dx*ix + dx
        ymin = dy*iy
        ymax = dy*iy + dy
        nmax = ss[ix, iy]
        realized_ss[ix, iy], grid_rows = sample_cell(individuals, xmin, xmax, ymin, ymax, nmax)
        sample_rows = np.concatenate([sample_rows, grid_rows])
    return(realized_ss, sample_rows)

def plot_samples(sample_parents, image_w, image_h, max_width, max_height):
    # creating new Image object
    w, h = image_w, image_h
    sampling = Image.new("1", (w, h))
    img2 = ImageDraw.Draw(sampling)

    # Plot sampling
    for ind_i, ind_row in sample_parents.iterrows():
        x, y = ind_row[['x','y']]
        img2.point((x*w/max_width, y*h/max_height), fill = "white")
    return(sampling)

def plot_intensity(sampling_intensity, image_w, image_h, max_int):
    w, h = image_w, image_h
    intensity = Image.new("L", (w, h))
    img3 = ImageDraw.Draw(intensity)

    # Width and height of grid cells
    x_cells = sampling_intensity.shape[0]
    y_cells = sampling_intensity.shape[1]
    dx = w/x_cells
    dy = h/y_cells
    # Maximum sampling intensity (will be 255 or all white)
    # Plot sampling intensity
    for ix, iy in np.ndindex(sampling_intensity.shape):
        xmin = dx*ix
        xmax = dx*ix + dx
        ymin = dy*iy
        ymax = dy*iy + dy
        # Sampling
        rel_intensity = round(255*sampling_intensity[ix, iy]/max_int)
        img3.rectangle([(xmin, ymin), (xmax, ymax)], fill = rel_intensity, outline=None, width=0)
    return(intensity)

def find_POPs_or_sibs(focal_parents, sample_array, type):
    # focal_parents: List of individual IDs of the focal parents. The focal parents do not need to be in the sample.
    # sample_array: Array of information about the whole sample
    # type: "PO" or "HS"
    #   "PO": find pairs with parent in focal_parent_array and offspring in sample_array
    #   "sib": find pairs of children in sample_array that share a parent in focal_parent_array
    # Returns:
    # An array of tuples of either [(parent, offspring), ...] or [(sibling 1, sibling 2), ...]
    relative_pairs = []
    for parent_id in focal_parents:
        #print("Parent", parent_id)
        children = sample_array[np.logical_or(sample_array.loc[:,'parent1'] == parent_id,
                                        sample_array.loc[:,'parent2'] == parent_id)]
        #print("Children")
        #print(children)
        if type == "PO":
        # Record POPs
            for j, child in children.iterrows():
                child_id = child['individual']
                relative_pairs.append((parent_id, child_id))
        if type == "HS":
            if len(children) == 2:
                sib_ids = children.loc[:,'individual'].values
                relative_pairs.append(tuple(np.sort(sib_ids)))
            if len(children) > 2:
                for pair in itertools.combinations(children.loc[:,'individual'], 2):
                    relative_pairs.append(tuple(np.sort(pair)))
    return(relative_pairs)

def draw_pairs(pairs, image, sample_array, max_width, max_height, color = "white"):
    # pairs: An array of tuples of individual ids
    # image: An image object to add lines to
    # sample_array: Array containing information about each individual
    for pair in pairs:
        index1 = np.where(sample_array.loc[:,'individual'] == pair[0])
        index2 = np.where(sample_array.loc[:,'individual'] == pair[1])
        row1 = sample_array.iloc[index1]
        row2 = sample_array.iloc[index2]
        x1, y1 = row1[['x', 'y']].values[0]
        x2, y2 = row2[['x', 'y']].values[0]
        image.line([(x1*w/max_width, y1*h/max_height), (x2*w/max_width, y2*h/max_height)], fill =color, width = 0)