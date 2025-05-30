# CKMRnn: Simulation-based spatially explicit close kin mark recapture
## Estimation of population size in African elephants

To replicate the results of the elephant analysis in our paper, run the following commands. More detailed descriptions of the scripts and data files are below.

1. Install and load required packages using [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).

```bash
cd elephants
mamba env create -f environment.yml
conda activate elephants
```

2. Process the empirical data.
```bash
Rscript processing.R
```

3. Generate training dataset.
```bash
cd simulations
snakemake -c8
```

4. Generate images from the empirical data
```bash
cd ../network/empirical_input_images
python empirical_images.py
```

5. Train network and generate parametric bootstrap replicates.
```bash
cd ..
snakemake -c8
```

6. Plot results
```
Rscript nn_results.R
```

### Empirical data

`elephants/empirical_data/Kibale_Species_Classification_Final.csv`: Information on all collected dung samples inculding individual ids, locations, and sampling dates. From Goodfellow et al. (2025).

`elephants/empirical_data/MLRelate_KNP.csv`: Kin relationships for all sampled individuals.

`elephants/empirical_data/Ug_Protected-areas2007/`: Shape file containing all protected areas in Uganda.

`elephants/processing.R`: Script to processes empirical data for use in the simulation-based CKMR analysis.

Output

`elephants/simulations/kibale.png`: Map of Kibale to use in the simulation.

`elephants/empirical_data/processed/`: Data frames with information about sampling, recaptures, and parent-offspring pairs in the empirical data.

### Simulation

`elephants/simulations/Snakefile`: Snakemake workflow to generate training dataset.

`elephants/simulations/config.yaml`: Config file for the snakemake workflow to set the number of training simulations to run and the parameter ranges for these simulations.

Output

`elephants/simulations/elephant_simulations/spatial_labels.csv`: Summary of the training dataset, including true population sizes, parameters used to run the simulations, and file locations for kin, recapture, and sampling intensity images.

### Network training, network testing, prediction, and parametric bootstrapping

`elephants/network/empirical_input_images/empirical_images.py`: Creates images of half sibling pairs, recaptures, and sampling intensity from the empirical data.

`elephants/network/Snakefile`: Snakemake workflow for network training and testing, predicton of population size for the empirical data, and parametric bootstrapping. Running the snakemake workflow requires access to a GPU.

`elephants/network/config.yaml`: Config file for the snakemake workflow specifying the set of simulations to use for training and the empirical images.

`elephants/network/nn_results.R`: Code for generating figures and computing 95% CI.

Output

`elephants/network/empirical_input_images`: Images of the empirical data for input to the neural network.

`elephants/network/network_output/training_hist.csv`: Training loss over the 20 training epochs.

`elephants/network/network_output/test_res.csv`: Estimated population sizes for the held out test simulations.

`elephants/network/network_output/pred.txt`: Estimated population size for the empirical data.

`elephants/network/network_output/bootstrap_reps.csv`: Estimated population sizes for the parametric bootstrap simulations.

### Reference
1. Goodfellow CK, Chusyd DE, Bird SR, Babaasa D, Chapman CA, Hickey JR, et al. Elephants inhabiting two forested sites in western Uganda exhibit contrasting patterns of species identity, density, and history of hybridization. 2025;:2025.05.13.653790.

## Simulation tests

To install and load required packages to replicate the simulation tests using [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html):

```bash
cd simulation_tests
mamba env create -f environment.yml
conda activate ckmr_sim_tests
```
Leslie matrix model for bearded seals: `leslie_model.R`

### Constant size simulations
`simulation_tests/bearded_seals/Snakefile`: Snakemake workflow to generate constant population size training dataset.

```bash
cd bearded_seals
snakemake --cores
```

Output

`simulation_test/bearded_seals/labels.csv`: Parameters, population sizes, and number of kin pairs for each training simulation.

`simulation_test/bearded_seals/bearded_seal_images/`: Kin and sampling intensity images for the training simulations

### Varying trend simulations
`simulation_tests/bearded_seals_trend/Snakefile`: Snakemake workflow to generate varying population trend training dataset.

```bash
cd bearded_seals_popsize_change
snakemake --cores
```

Output

`simulation_tests/bearded_seals_trend/labels.csv`: Parameters, population sizes, and number of kin pairs for each training simulation.

`simulation_tests/bearded_seals_trend/bearded_seal_images/`: Kin and sampling intensity images for the training simulations

### Network training and testing
#### Accuracy with biased sampling

`simulation_tests/constant_size_network.ipynb`: 

Output

`simulation_tests/bearded_seals_nn_results/constant_size_network.csv`: Performance of the trained network on held out test simulations. 

#### The effect of misspecified population trend

`simulation_tests/bearded_seals_trend/misspecified_trend.ipynb`

Output

`simulation_tests/bearded_seals_nn_results/misspecified_trend.csv`: Performance of network trained on constant size simulations when tested on simulations with varying population trends.

#### Training for robustness to population trend

`simulation_tests/bearded_seals_trend/varying_trend_network.ipynb`

Output

`simulation_tests/bearded_seals_nn_results/varying_trend_network.csv`: Performance of network trained on varying population trend simulations when tested on a held out set of simulations with varying population trends.

### Visualizing results

`simulation_tests/bearded_seals_nn_results/nn_plots.R`