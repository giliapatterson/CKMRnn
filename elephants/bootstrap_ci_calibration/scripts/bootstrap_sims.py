import os
import numpy as np
import pandas as pd
import subprocess

simulation_folder = snakemake.config['simulation_folder']
workflow_path = f"{snakemake.config['simulation_folder']}/Snakefile"
configfile = f"{snakemake.config['simulation_folder']}/config.yaml"
cores = 20

# Read in predicted population sizes
predictions = pd.read_csv(snakemake.input.pred)

# Seed
rng = np.random.default_rng(snakemake.params.seed)

N_pred_unrounded = predictions['pred']
predictions['N_pred'] = round(N_pred_unrounded, 0).astype(int)
print(predictions)
# Find N_pred that appear more than once, and create unique output folders for each
num_copies = predictions['N_pred'].value_counts()
# Suffix to add to folder names
predictions['num_copies'] = predictions['N_pred'].map(num_copies)
suffix = ['' if predictions['num_copies'][i]==1 else f"_{i+1}" for i in range(predictions.shape[0])]

# Run bootstrap sims for each prediction
for i, row in predictions.iterrows():
    seed = rng.integers(1, 1e6)
    output_folder = f"{snakemake.config['simulation_folder']}/bootstrap_sims_N{row['N_pred'].astype(int)}{suffix[i]}/"
    command = f"snakemake --cores {cores} -s {workflow_path} --configfile {configfile} -d {simulation_folder} --config folder={output_folder} seed={seed} minN={row['N_pred']} maxN={row['N_pred']} nreps={snakemake.params.nreps}"
    # Run the Snakemake workflow
    retcode = subprocess.call(command, shell=True)
    if retcode != 0:
        print(f"Simulation for N={row['N_pred']} failed, skipping")
        predictions.loc[i, 'bootstrap_folder'] = 'none'
        continue
    predictions.loc[i, 'bootstrap_folder'] = output_folder
predictions.to_csv(snakemake.output.predictions_with_bootstrap, index = False)
