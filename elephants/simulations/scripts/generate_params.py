import pickle
import numpy as np
import pandas as pd
rng = np.random.default_rng(snakemake.params.seed)

num_combs = len(snakemake.params.reps)
# Draw N from the prior distribution
N = rng.uniform(snakemake.params.minN, snakemake.params.maxN, size = num_combs).astype(int)
# Draw dispersal and mate choice parameters from uniform distributions
DISPERSAL_SIGMA = rng.uniform(snakemake.params.min_DS, snakemake.params.max_DS, size = num_combs)
SIGMA = rng.uniform(snakemake.params.min_MC, snakemake.params.max_MC, size = num_combs)

# Draw seeds for SLiM simulations
SEED = rng.integers(low = 1, high = 3000, size = num_combs).astype(int)

all_combinations = {"N": N,
                    "n": snakemake.params.n,
                    "rep": snakemake.params.reps,
                    "DISPERSAL_SIGMA": DISPERSAL_SIGMA,
                    "SIGMA": SIGMA,
                    "SEED": SEED,
                    "max_DS": snakemake.params.max_DS,
                    "min_DS": snakemake.params.min_DS,
                    "max_MC": snakemake.params.max_MC,
                    "min_MC": snakemake.params.min_MC,
                    "FAGEREPRO": snakemake.params.FAGEREPRO,
                    "MAGEREPRO": snakemake.params.MAGEREPRO,
                    "SURVIVAL_POWER": snakemake.params.SURVIVAL_POWER,
                    "SAMPLE_SIGMA": snakemake.params.SAMPLE_SIGMA,
                    "PSAMPLE": snakemake.params.PSAMPLE,
                    "MAXWIDTH": snakemake.params.MAXWIDTH,
                    "MAXHEIGHT": snakemake.params.MAXHEIGHT,
                    "DW": snakemake.params.DW,
                    "DH": snakemake.params.DH,
                    "pxperkm": snakemake.params.pxperkm
                    }
        

df = pd.DataFrame(all_combinations)      

print("data frame generated")
# chop up the df
out_files = snakemake.output.df
breakpoints = np.linspace(0, df.shape[0], len(out_files) + 1).astype(int)
print(breakpoints)
for output, start, end in zip(out_files, breakpoints[:-1], breakpoints[1:]):
    subset_df = df[start:end]
    pickle.dump(subset_df, open(output, 'wb'))