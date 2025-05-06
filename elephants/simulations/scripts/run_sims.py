import os
import numpy as np
import pandas as pd
import pickle

# Load the data frame
params_df = pickle.load(open(snakemake.input.df, 'rb'))
# Add columns for file names
params_df['spatial_sample'] = ''
params_df['random_sample'] = ''
params_df['spatial_sample_locations'] = ''
params_df['popsize_file'] = ''

folder = f"{snakemake.params.folder}slim_output_{snakemake.params.BATCH_ID}"
os.mkdir(folder)

# Loop through rows and run simulations
for i, row in params_df.iterrows():

    # Generate output file names
    spatial_sample = f"{folder}/N{row['N']}_{i}_spatial_sample.csv"
    random_sample = f"{folder}/N{row['N']}_{i}_random_sample.csv"
    spatial_sample_locations = f"{folder}/N{row['N']}_{i}_spatial_sample_locations.csv"
    popsize_file = f"{folder}/N{row['N']}_{i}_popsizes.csv"

    # Run SLiM
    command = f"slim -d RANDOM_SEED={row['SEED']} -d 'MAP_FILE=\"{snakemake.input.kibale_image}\"' -d 'SURVIVAL_FILE=\"{snakemake.input.survival_file}\"' -d 'SAMPLE_TIMES_FILE=\"{snakemake.input.sample_times}\"' -d 'POTENTIAL_LOCS_FILE=\"{snakemake.input.sampling_locations}\"' -d N0={row['N']} -d SAMPLE_SIZE={row['n']} -d DISPERSAL_SIGMA={row['DISPERSAL_SIGMA']} -d SIGMA={row['SIGMA']} -d FAGEREPRO={row['FAGEREPRO']} -d MAGEREPRO={row['MAGEREPRO']} -d SURVIVAL_POWER={row['SURVIVAL_POWER']} -d SAMPLE_SIGMA={row['SAMPLE_SIGMA']} -d PSAMPLE={row['PSAMPLE']} -d MAX_WIDTH={row['MAXWIDTH']} -d MAX_HEIGHT={row['MAXHEIGHT']} -d 'SPATIAL_SAMPLE=\"{spatial_sample}\"' -d 'SPATIAL_SAMPLE_LOCS=\"{spatial_sample_locations}\"' -d 'RANDOM_SAMPLE=\"{random_sample}\"' -d 'POPSIZE_FILE=\"{popsize_file}\"' elephants.slim"
    print(command)
    os.system(command)

    # Add output file names to the dataframe
    params_df.loc[i,'spatial_sample'] = spatial_sample
    params_df.loc[i,'random_sample'] = random_sample
    params_df.loc[i,'spatial_sample_locations'] = spatial_sample_locations
    params_df.loc[i,'popsize_file'] = popsize_file 

pickle.dump(params_df, open(snakemake.output.df, 'wb'))