import os
import numpy as np
import pandas as pd
import os

simulation_folder = snakemake.config['simulation_folder']
workflow_path = f"{snakemake.config['simulation_folder']}/Snakefile"
configfile = f"{snakemake.config['simulation_folder']}/config.yaml"
cores = 8

# Folder to store simulations
output_folder = f"{snakemake.config['simulation_folder']}/test_sims_for_CI_calibration/"

command = f"snakemake --cores {cores} -s {workflow_path} --configfile {configfile} -d {simulation_folder} --config folder={output_folder} seed={snakemake.params.seed} nreps={snakemake.params.nreps}"

# Run the Snakemake workflow
os.system(command)
with open(snakemake.output.output_folder, "w") as file:
    file.write(output_folder)
