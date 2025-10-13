import numpy as np
import pandas as pd

import sample_and_estimate as se

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--labels', type=str)
parser.add_argument('--base_folder', type=str, default = "../bearded_seals")
parser.add_argument('--noptsteps', type=int, default=50)
parser.add_argument('--output_file', type=str, default="bearded_seal_estimates.csv")
args = parser.parse_args()

labels = pd.read_csv(args.labels)
base_folder = args.base_folder
labels = labels.iloc[0:9,]
noptsteps = args.noptsteps

for i, row in labels.iterrows():
    R0 = row.loc['R0']
    rep = int(row['rep'])
    bias = row['bias']
    # If row contains surv_mult, read it in
    if 'surv_mult' in row.index:
        surv_mult = row.loc['surv_mult']
        parents_path = f'{base_folder}/bearded_seal_parents/parents_{R0}_{surv_mult}_{rep}.csv'
    else:
        parents_path = f'{base_folder}/bearded_seal_parents/parents_{R0}_{rep}.csv'
    print(f"Estimating for simulation {i + 1} of {labels.shape[0]}")
    estR0, estN = se.sample_and_estimate(parents_path, bias, noptsteps)
    labels.loc[i, 'estR0'] = estR0
    labels.loc[i, 'estN'] = estN

labels.to_csv(args.output_file)
