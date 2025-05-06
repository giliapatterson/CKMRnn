import numpy as np
from scipy import stats
import sys
import itertools
import math
import pandas as pd
from PIL import Image, ImageDraw
import runpy
import sampling_and_kin_functions as skf

# input
pops = pd.read_csv("../../empirical_data/processed/empirical_pops_df.csv")
#sibs = pd.read_csv("../../empirical_data/processed/empirical_sibs_df.csv")
recaps = pd.read_csv("../../empirical_data/processed/empirical_recaps_df.csv")
spatial_sample_intensity = pd.read_csv("../../empirical_data/processed/empirical_intensity_df.csv")

# Parameters from the config file
max_width = 37.76
max_height = 55.85
DW = 1.021
DH =  1.015
w = int(max_width * 15)
h = int(max_height * 15)

# Plot sampling intensity
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
intensity.save("empirical_intensity.png")

# Plot pops, sibs, and recaptures
pops_spaghetti = Image.new("1", (w, h))
img1 = ImageDraw.Draw(pops_spaghetti)
skf.draw_pairs(pops, img1, max_width, max_height, w, h)

sibs_spaghetti = Image.new("1", (w, h))
img2 = ImageDraw.Draw(sibs_spaghetti)
skf.draw_pairs(sibs, img2, max_width, max_height, w, h)

recaptures_spaghetti = Image.new("1", (w, h))
img3 = ImageDraw.Draw(recaptures_spaghetti)
skf.draw_pairs(recaps, img3, max_width, max_height, w, h)

pops_spaghetti.save("empirical_pops.png")
#sibs_spaghetti.save("empirical_sibs.png")
recaptures_spaghetti.save("empirical_recaps.png")