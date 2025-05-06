import pickle
import pandas as pd

dfs = [pickle.load(open(file, 'rb')) for file in snakemake.input.df]
df = pd.concat(dfs)
df.to_csv(snakemake.output.labels, index=False)