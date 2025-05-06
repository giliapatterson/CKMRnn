import numpy as np

import torch
from torch.utils.data import Dataset
import pandas as pd
from torch.utils.data import DataLoader, random_split
from torch import nn
from torchvision.io import read_image

ptrain = snakemake.params.ptrain

torch.manual_seed(100)
rng = np.random.default_rng(seed = 200)

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using {device} device", flush=True)

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
    
spatial_labels = pd.read_csv(snakemake.input.labels)
sim_folder = snakemake.params.folder

kin_dataset = KinDataset(spatial_labels, sim_folder)

n_total = len(kin_dataset)
n_train = round(n_total*ptrain)
n_test = n_total-n_train

train_kin, test_kin = random_split(kin_dataset, [n_train, n_test], generator=torch.Generator().manual_seed(42))

batch_size = 10
train_kin_dl = DataLoader(train_kin, batch_size, shuffle = True)
test_kin_dl = DataLoader(test_kin, batch_size, shuffle = True)

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

model = Net()
model = model.to(device) 

loss_fn = nn.MSELoss(reduction='mean')

optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)

def train(model, num_epochs, train_dl):
    loss_hist_train = [0] * num_epochs
    for epoch in range(num_epochs):
        model.train()
        for x_batch, y_batch, ds, mc in train_dl:
            x_batch = x_batch.to(device) 
            y_batch = y_batch.to(device) 
            pred = model(x_batch)
            loss = loss_fn(pred, y_batch)
            loss.backward()
            optimizer.step()
            optimizer.zero_grad()
            loss_hist_train[epoch] += loss.item()*y_batch.size(0)

        loss_hist_train[epoch] /= len(train_dl.dataset)
                
        print(f'Epoch {epoch+1} loss: {loss_hist_train[epoch]:.4f}')
    return loss_hist_train

torch.manual_seed(1)
num_epochs = 20
hist = train(model, num_epochs, train_kin_dl)

#Save model
PATH = snakemake.output.network
torch.save(model.state_dict(), PATH)

train_x = np.arange(len(hist))
hist_df = pd.DataFrame({'train_x': train_x, 'train_loss': hist})
hist_df.to_csv(snakemake.output.hist, index = False)

test_kin_small, extra_test_kin = random_split(test_kin, [100, n_test-100], generator=torch.Generator().manual_seed(400))
torch.cuda.synchronize()
model_cpu = model.cpu()

test_data = test_kin
test_truth = np.empty(len(test_data))
test_pred = np.empty(len(test_data))
test_ds = np.empty(len(test_data))
test_mc = np.empty(len(test_data))

model.eval()
with torch.no_grad():
    for i, (test_input, test_output, ds, mc) in enumerate(test_data):
        test_pred[i] = model(test_input.unsqueeze(0))[0][0]
        test_truth[i] = test_output[0]
        test_ds[i] = ds
        test_mc[i] = mc

# Write test results to file
results = {'N': test_truth, 'pred': test_pred, 'ds': test_ds, 'mc': test_mc}
df = pd.DataFrame(data=results)
df.to_csv(snakemake.output.test_res, index = False)

# Get estimates from empirical data
intensity = read_image(snakemake.input.empirical_intensity)
pops = read_image(snakemake.input.empirical_pops)
recaps = read_image(snakemake.input.empirical_recaps)
empirical_input = torch.cat((pops, recaps, intensity), 0).float()
model.eval()
N_pred = model(empirical_input.unsqueeze(0))[0][0].item()

with open(snakemake.output.pred, "w") as file:
    file.write(str(N_pred))