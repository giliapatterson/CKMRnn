import numpy as np

import torch
from torch.utils.data import Dataset
import pandas as pd
from torch.utils.data import DataLoader, random_split
from torch import nn
from torchvision.io import read_image

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using {device} device")

class KinDataset(Dataset):
    def __init__(self,labels_dataframe,sim_folder):
        self.labels = labels_dataframe
        self.sim_folder = sim_folder
    def __len__(self):
        return len(self.labels)
    def __getitem__(self, idx):
        row = self.labels.iloc[idx,]
        spaghetti_pops = read_image(f"{self.sim_folder}/{row.loc['spaghetti_pops']}")
        spaghetti_recaps = read_image(f"{self.sim_folder}/{row.loc['spaghetti_recaptures']}")
        intensity_map = read_image(f"{self.sim_folder}/{row.loc['intensity_image']}")
        
        input_tensor = torch.cat((spaghetti_pops, spaghetti_recaps, intensity_map), 0).float()
        label = torch.tensor([self.labels['N_final'].iloc[idx]]).float()
        ds = row.loc['DISPERSAL_SIGMA']
        mc = row.loc['SIGMA']
        return input_tensor, label, ds, mc

class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.conv1 = nn.Conv2d(in_channels=3, out_channels=32, kernel_size=6, padding=3)
        self.conv2 = nn.Conv2d(in_channels=32, out_channels=64, kernel_size=6, padding=3)
        self.fc1 = nn.Linear(1908480, 1024)
        self.fc2 = nn.Linear(1024, 1)
        self.dropout = nn.Dropout(p=0.2)
        self.pool = nn.MaxPool2d(kernel_size=2)
        self.relu = nn.ReLU()
        self.flatten = nn.Flatten() 
    def forward(self, x):
        x = self.pool(self.relu(self.conv1(x)))
        x = self.pool(self.relu(self.conv2(x)))
        x = self.flatten(x)
        x = self.dropout(self.relu(self.fc1(x)))
        x = self.fc2(x)
        return x

# Load model
PATH = snakemake.input.network
model = Net()
model.load_state_dict(torch.load(PATH))
model.eval()

# Predict on each bootstrap replicate
torch.cuda.synchronize()
model_cpu = model.cpu()

# Read in data to predict on
with open(snakemake.input.sim_folder, 'r') as file:
    sim_folder = file.readline().strip()
spatial_labels = pd.read_csv(f"{sim_folder}spatial_labels.csv")
bootstrap_dataset = KinDataset(spatial_labels, '.')

test_kin = bootstrap_dataset
test_pred = np.empty(len(test_kin))

with torch.no_grad():
    for i, (test_input, test_output, ds, mc) in enumerate(test_kin):
        test_pred[i] = model(test_input.unsqueeze(0))[0][0]

# Write bootstrap results to file
results = {'pred': test_pred, 'estimate': spatial_labels.loc[:,'N']}
df = pd.DataFrame(data=results)
df.to_csv(snakemake.output.bootstrap_reps, index = False)