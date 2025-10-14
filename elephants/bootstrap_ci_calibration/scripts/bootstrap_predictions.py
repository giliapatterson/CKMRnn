import os
import numpy as np
import pandas as pd
import subprocess

# Read in test simulation information, including the folders where bootstrap simulations are located
test_sims = pd.read_csv(snakemake.input.pred)

# Folder to store bootstrap predictions
bootstrap_pred_folder = snakemake.params.bootstrap_pred_folder
if not os.path.exists(bootstrap_pred_folder):
    os.makedirs(bootstrap_pred_folder)
# Predict on bootstrap sims for each simulation
for i, row in test_sims.iterrows():
    # Predict on bootstrap sims
    N_pred = round(row['pred'])
    output_file = f"{bootstrap_pred_folder}/predictions_N{N_pred}.csv"
    if os.path.exists(output_file):
        print("Predictions already exist, skipping")
        continue
    command = f"python scripts/predict.py --sim_folder {row['bootstrap_folder']} --network {snakemake.input.network} --predictions_file {output_file}"
    retcode = subprocess.call(command, shell=True)
    # Check if the simulation finished, and if not, don't record any values in the dataframe
    if retcode == 0:
        test_sims.loc[i, 'bootstrap_predictions_file'] = output_file

test_sims.to_csv(snakemake.output.predictions_with_bootstrap_reps, index = False)
