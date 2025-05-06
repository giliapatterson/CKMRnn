# CKMRnn: Simulation-based spatially explicit close kin mark recapture
## Simulation tests

To install and load required packages to replicate the simulation tests:

```bash
cd simulation_tests
conda env create -f environment.yml
conda activate ckmr_sim_tests
```
Leslie matrix model to calculate parameters for constant population size: `leslie_model.R`

### Constant size simulations
`simulation_tests/bearded_seals`

### Varying trend simulations
`simulation_tests/bearded_seals_popsize_change`

### Network training and testing
#### Accuracy with biased sampling
`simulation_tests/bearded_seal_network_sibs.ipynb`
"bearded_seals_nn_results/model_pops_sibs.csv"

#### The effect of misspecified population trend
simulation_tests/bearded_seals_popsize_change/test_sibling_network.ipynb
bearded_seals_nn_results/model_pops_sibs_changing size.csv

#### Training for robustness to population trend
simulation_tests/bearded_seals_popsize_change/bearded_seal_network_sibs_popsize_decline.ipynb
/bearded_seals_nn_results/model_pops_sibs_popsize_change.csv

## Estimation of population size in African elephants

To install and load required packages to replicate the African elephant CKMRnn analysis:

```bash
cd elephants
conda env create -f environment.yml
conda activate elephants
```

### Empirical data

`elephants/empirical_data/Kibale_Species_Classification_Final.csv`: Information on all collected dung samples inculding individual ids, locations, and sampling dates. From Goodfellow et al. (2025).

`elephants/empirical_data/MLRelate_KNP.csv`: Kin relationships for all sampled individuals.

`elephants/empirical_data/Ug_Protected-areas2007/`: Shape file containing all protected areas in Uganda.

`elephants/processing.R`: Script to processes empirical data for use in the simulation-based CKMR analysis.

### Simulation

`elephants/simulations/elephant_simulations/spatial_labels.csv`: Summary of the simulations used for training and testing the neural network, including true population size and parameters used to run the simulation.

Training data for the neural network was generated using a snakemake workflow, with the number of simulations to run and the parameter ranges set in `config.yaml`. To run the workflow:

```bash
cd simulations
snakemake -c8
```

### Network training, network testing, prediction, and parametric bootstrapping

`network/empirical_input_images`: Images of the empirical data for input to the neural network. Generated using the script `empirical_images.py`.

`network/network_output/training_hist.csv`: Training loss over the 20 training epochs.

`network/network_output/test_res.csv`: Estimated population sizes for the held out test simulations.

`network/network_output/pred.txt`: Estimated population size for the empirical Kibale elephants (450).

`network/network_output/bootstrap_reps.csv`: Estimated population sizes for the parametric bootstrap simulations.

`network/nn_results.R`: Figures.

Network training and testing, predicton of population size for the empirical data, and parametric bootstrapping were all performed using a Snakemake workflow. Running the workflow requires access to a GPU. To run the workflow:

```bash
cd network
snakemake -c8
```