import numpy as np
import pandas as pd

# Read in test simulation information, including the folders where bootstrap replicates are located
test_sims = pd.read_csv(snakemake.input.pred)

# Itterate through test simulations and compute rank of the true value
for i, row in test_sims.iterrows():
    bootstrap_predictions_file = row['bootstrap_predictions_file']
    print(f"Reading in {bootstrap_predictions_file}")
    bootstrap_predictions = pd.read_csv(bootstrap_predictions_file)
    true_value = row['truth']
    rank = sum(bootstrap_predictions['pred'] < true_value) + sum(bootstrap_predictions['pred'] == true_value)/2
    test_sims.loc[i, 'rank'] = rank
    test_sims.loc[i, 'n_bootstrap_reps'] = bootstrap_predictions.shape[0]
    test_sims.loc[i, 'rank_normalized'] = rank/bootstrap_predictions.shape[0]
    # Is the true value in the 95% CI?
    lower = np.percentile(bootstrap_predictions['pred'], 2.5)
    upper = np.percentile(bootstrap_predictions['pred'], 97.5)
    test_sims.loc[i, 'in_95_CI'] = (true_value >= lower) & (true_value <= upper)
    # Is the true value in the 50% CI?
    lower = np.percentile(bootstrap_predictions['pred'], 25)
    upper = np.percentile(bootstrap_predictions['pred'], 75)
    test_sims.loc[i, 'in_50_CI'] = (true_value >= lower) & (true_value <= upper)
test_sims.to_csv(snakemake.output.predictions_with_ranks, index = False)

