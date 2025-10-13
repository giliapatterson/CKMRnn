import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--paths', nargs='+', type=str)
parser.add_argument('--output_file', type=str, default="all_estimates.csv")
args = parser.parse_args()

dfs = [pd.read_csv(file) for file in args.paths]
df = pd.concat(dfs)
df.to_csv(args.output_file, index=False)
