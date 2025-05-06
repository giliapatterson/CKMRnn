
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
# Add column for file names
params_df['intensity_image'] = ''

# Make folder for output files
folder = f"{snakemake.params.folder}intensity_{snakemake.params.BATCH_ID}"
os.mkdir(folder)

# Loop through rows and plot intensity images
for i, row in params_df.iterrows():
    intensity_in = row['spatial_sample_locations']
    image_out = f"{folder}/N{row['N']}_{i}_intensity.png"
    dataframe_out = f"{folder}/N{row['N']}_{i}_intensity_df.csv"
    params_df.loc[i,'intensity_image'] = image_out

    max_width = row['MAXWIDTH']
    max_height = row['MAXHEIGHT']
    DW = row['DW']
    DH = row['DH']

    spatial_sample_intensity = pd.read_csv(intensity_in)
    spatial_sample_intensity = spatial_sample_intensity.groupby(['x','y'])['nsampled'].sum().reset_index()
    spatial_sample_intensity.to_csv(dataframe_out, index = False)

    w = int(max_width * row['pxperkm'])
    h = int(max_height * row['pxperkm'])

    intensity = Image.new("L", (w, h))
    img1 = ImageDraw.Draw(intensity)
    xlist = spatial_sample_intensity.loc[:,'x'].values
    ylist = max_height - spatial_sample_intensity.loc[:,'y'].values
    fills = spatial_sample_intensity.loc[:,'nsampled'].values.astype(int)

    maxfill = max(fills)
    for i, fill in enumerate(fills):
        x0 = (xlist[i] - DW/2)*w/max_width
        y0 = (ylist[i] - DH/2)*h/max_height
        x1 = (xlist[i] + DW/2)*w/max_width
        y1 = (ylist[i] + DH/2)*h/max_height
        fill_n = int(255*fill/maxfill)
        img1.rectangle([x0, y0, x1, y1], fill = fill_n)
    intensity.save(image_out)

# Write dataframe to file
pickle.dump(params_df, open(snakemake.output.df, 'wb'))