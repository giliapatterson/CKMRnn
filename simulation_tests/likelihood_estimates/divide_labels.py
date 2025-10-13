import pandas as pd
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--paths', nargs='+', type=str)
parser.add_argument('--input_file', type=str, default = "../bearded_seals_trend/labels.csv", help='Path to the input labels CSV file')
args = parser.parse_args()

paths = args.paths
ndivide = len(paths)
labels = pd.read_csv(args.input_file)
#print(paths)

# Divide labels into ndivide parts
labels_split = np.array_split(labels, ndivide)

# Save each part to a separate CSV files
for i, part in enumerate(labels_split):
    #print(paths[i])
    part.to_csv(paths[i], index=False)
