import os
import numpy as np
import pandas as pd
import os

simulation_folder = snakemake.config['simulation_folder']
workflow_path = f"{snakemake.config['simulation_folder']}/Snakefile"
configfile = f"{snakemake.config['simulation_folder']}/config.yaml"
cores = 8

# Read in predicted population size
with open(snakemake.input.pred, 'r') as file:
    line = file.readline().strip()
N_pred = round(float(line))

# Folder to store simulations
output_folder = f"{snakemake.config['simulation_folder']}/bootstrap_sims_N{N_pred}/"

command = f"snakemake -s {workflow_path} --configfile {configfile} -d {simulation_folder} --config folder={output_folder} seed={snakemake.params.seed} minN={N_pred} maxN={N_pred} nreps={snakemake.params.nreps}"

# Run the Snakemake workflow
os.system(command)
with open(snakemake.output.output_folder, "w") as file:
    file.write(output_folder)
