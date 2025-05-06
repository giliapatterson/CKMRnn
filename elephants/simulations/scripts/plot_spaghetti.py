import numpy as np
from scipy import stats
import sys
import itertools
import math
import pandas as pd
from PIL import Image, ImageDraw
import sampling_and_kin_functions as skf
import os
import pickle

# Load the data frame
params_df = pickle.load(open(snakemake.input.df, 'rb'))
# Add columns for file names
params_df['spaghetti_pops'] = ''
params_df['spaghetti_sibs'] = ''
params_df['spaghetti_recaptures'] = ''
# Add columns for population size and number of pairs
params_df['N_final'] = ''
params_df['N_avg'] = ''
params_df['npops'] = ''
params_df['nsibs'] = ''
params_df['nrecaps'] = ''

# Folders for output files
folder = f"{snakemake.params.folder}spaghetti_images_{snakemake.params.BATCH_ID}"
os.mkdir(folder)

# Loop through rows and plot intensity images
for i, row in params_df.iterrows():
    # Input files
    parents_file = row['spatial_sample']
    popsize_file = row['popsize_file']

    # Output files
    spaghetti_pops_out = f"{folder}/N{row['N']}_{i}_spaghetti_pops.png"
    spaghetti_sibs_out = f"{folder}/N{row['N']}_{i}_spaghetti_sibs.png"
    spaghetti_recaptures_out = f"{folder}/N{row['N']}_{i}_spaghetti_recaptures.png"
    pops_df = f"{folder}/N{row['N']}_{i}_pops_df.csv"
    sibs_df = f"{folder}/N{row['N']}_{i}_sibs_df.csv"
    recaptures_df = f"{folder}/N{row['N']}_{i}_recaptures_df.csv"

    # Add columns for file names
    params_df.loc[i,'spaghetti_pops'] = spaghetti_pops_out
    params_df.loc[i,'spaghetti_sibs'] = spaghetti_sibs_out
    params_df.loc[i, 'spaghetti_recaptures'] = spaghetti_recaptures_out

    # Sampled individuals
    sample_parents = pd.read_csv(parents_file)

    # Record POPs, recaptures, and half-sibling pairs
    pops = skf.find_pops(sample_parents)
    pops.to_csv(pops_df, index = False)
    npops = len(pops)

    hs = skf.find_sibs(sample_parents)
    hs.to_csv(sibs_df, index = False)
    nsibs = len(hs)

    recaptures = skf.find_recaptures(sample_parents)
    recaptures.to_csv(recaptures_df, index = False)
    nrecaps = len(recaptures)

    # Plot maternal POPs, paternal POPs, and half-sibling pairs
    max_width = row['MAXWIDTH']
    max_height = row['MAXHEIGHT']
    w = int(max_width * row['pxperkm'])
    h = int(max_height * row['pxperkm'])
    pops_spaghetti = Image.new("1", (w, h))
    img1 = ImageDraw.Draw(pops_spaghetti)
    skf.draw_pairs(pops, img1, max_width, max_height, w, h)

    sibs_spaghetti = Image.new("1", (w, h))
    img2 = ImageDraw.Draw(sibs_spaghetti)
    skf.draw_pairs(hs, img2, max_width, max_height, w, h)

    recaptures_spaghetti = Image.new("1", (w, h))
    img3 = ImageDraw.Draw(recaptures_spaghetti)
    skf.draw_pairs(recaptures, img3, max_width, max_height, w, h)

    pops_spaghetti.save(spaghetti_pops_out)
    sibs_spaghetti.save(spaghetti_sibs_out)
    recaptures_spaghetti.save(spaghetti_recaptures_out)


    # Population sizes
    popsize = pd.read_csv(popsize_file)
    # Calculate average population size
    N_avg = np.mean(popsize.loc[:,'N'])
    N_final = popsize.loc[:, 'N'].values[-1]
    
    ### Writing metadata file ###
    params_df.loc[i,'N_final'] = N_final
    params_df.loc[i,'N_avg'] = N_avg
    params_df.loc[i,'npops'] = npops
    params_df.loc[i,'nsibs'] = nsibs
    params_df.loc[i,'nrecaps'] = nrecaps

# Write dataframe to file
pickle.dump(params_df, open(snakemake.output.df, 'wb'))
