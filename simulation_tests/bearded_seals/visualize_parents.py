import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
import sys
import itertools
import math
import pandas as pd
from PIL import Image, ImageDraw
rng = np.random.default_rng()

# R0, rep, bias
R0 = sys.argv[1]
rep = sys.argv[2]
bias = float(sys.argv[3])
max_bias = float(sys.argv[4])
max_n = float(sys.argv[5])
# Get the parents file path from the command line
parents_file = sys.argv[6]
popsize_file = sys.argv[7]
# Output files
spaghetti_pops_out = sys.argv[8]
spaghetti_sibs_out = sys.argv[9]
sampling_out = sys.argv[10]
intensity_out = sys.argv[11]
metadata_out = sys.argv[12]

# Parents of all dead individuals
parents = pd.read_csv(parents_file)

# Population sizes
popsize = pd.read_csv(popsize_file)
# Calculate average population size
N = np.mean(popsize.loc[:,'N'])


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

def plot_intensity(sampling_intensity, image_w, image_h, max_int = max_bias):
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

# Define sampling intensity grid with one side [bias] times as likely to be sampled as the other
sampling_intensity =  np.repeat([np.linspace(1, bias, 10)], 10, axis = 0)

# Sample according to the grid and return realized sampling size and sampled individuals
ss, sample_rows = sample_grid(parents, sampling_intensity, max_n, 10, 10)
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
maternal_pops = find_POPs_or_sibs(mother_list, sample_parents, "PO")
paternal_pops = find_POPs_or_sibs(father_list, sample_parents, "PO")
maternal_sibs = find_POPs_or_sibs(maternal_sib_parents, sample_parents, "HS")
paternal_sibs = find_POPs_or_sibs(paternal_sib_parents, sample_parents, "HS")

# Full sibling pairs appear in both maternal and paternal arrays
if np.any([x in paternal_sibs for x in maternal_sibs]):
    full_sibs = maternal_sibs[np.where([x in paternal_sibs for x in maternal_sibs])[0][0]]
else:
    full_sibs = []

# Plot maternal POPs, paternal POPs, and half-sibling pairs
w, h = 500, 500
max_width = 10
max_height = 10

pops_spaghetti = Image.new("1", (w, h))
img1 = ImageDraw.Draw(pops_spaghetti)
draw_pairs(maternal_pops, img1, sample_parents, max_width, max_height)
draw_pairs(paternal_pops, img1, sample_parents, max_width, max_height)

sibs_spaghetti = Image.new("1", (w, h))
img2 = ImageDraw.Draw(sibs_spaghetti)
draw_pairs(maternal_sibs, img2, sample_parents, max_width, max_height)
draw_pairs(paternal_sibs, img2, sample_parents, max_width, max_height)

# Number of parent offspring and HS pairs
npops = len(maternal_pops) + len(paternal_pops)
nsibs =  len(maternal_sibs) + len(paternal_sibs)

# Plot locations of samples, and sampling intensity
sampling = plot_samples(sample_parents, 500, 500, 10, 10)
intensity = plot_intensity(sampling_intensity, 500, 500)

pops_spaghetti.save(spaghetti_pops_out)
sibs_spaghetti.save(spaghetti_sibs_out)

sampling.save(sampling_out)
intensity.save(intensity_out)

n = np.sum(ss)

with open(metadata_out, "w") as f:
    f.write(str(R0) + "," + str(rep) + ","+ str(bias) + ","+ str(N) + "," + str(n) + "," + str(npops) + "," + str(nsibs) + "\n")