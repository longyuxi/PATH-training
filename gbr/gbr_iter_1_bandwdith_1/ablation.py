# %%
import numpy as np

# %%
from pathlib import Path
import pickle

num_datapoints_per_file = 500

# Loading the files back in
observations = []
binding_affinities = []

NPYS_FOLDER = Path('/usr/project/dlab/Users/jaden/gbr-tnet/persistence/persistence_images')

for file_index in range(len(list(NPYS_FOLDER.glob('binding_affinities_*.npy')))):
    with open(NPYS_FOLDER / f'observations_{file_index}.npy', 'rb') as f:
        observations.append(np.load(f))

    with open(NPYS_FOLDER / f'binding_affinities_{file_index}.npy', 'rb') as f:
        binding_affinities.append(np.load(f))


# %%
observations = np.concatenate(observations)
binding_affinities = np.concatenate(binding_affinities)

print(observations.shape)
print(binding_affinities.shape)


# %%
observations = observations.reshape(observations.shape[0], -1)
observation_indices = np.arange(observations.shape[1])

# %%
highly_selected_observations_indices = [16829, 16930, 17031, 17818, 17819, 17917, 17918, 17919, 17920, 18019, 18020, 18021, 18115, 18116, 18117, 18120, 18121, 18214, 18215, 18216, 18217, 18218, 18313, 18314, 18315, 18316, 18317, 18318, 18413, 18415, 18416, 18417, 18418, 18514, 18516, 18517, 18610, 18615, 18710, 18711, 18811, 18812, 18908, 18912, 18913, 19110, 287524, 287620, 287719, 287720, 287721, 287818, 287821, 287822, 287919, 287921, 287922, 287923, 288020, 288022, 288121, 288219, 288314, 288415, 288516, 31500, 31700, 31701, 31702, 61301, 91000, 91001, 91002, 91802, 92100, 92101, 92102]

highly_selected_observations = observations[:, highly_selected_observations_indices]
highly_selected_observations_indices = observation_indices[highly_selected_observations_indices]

# %%
print(highly_selected_observations.shape, highly_selected_observations_indices.shape, sep='\n')

# %%
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
import plotly.express as px
from sklearn.metrics import r2_score

# %%
stats = []

from tqdm import tqdm
from multiprocessing import Pool, cpu_count

Path('additional_results').mkdir(exist_ok=True)

def trial(i):
    X_train, X_test, y_train, y_test = train_test_split(highly_selected_observations, binding_affinities, test_size=0.1, random_state=i)
    regr = GradientBoostingRegressor(random_state=2*i)
    regr.fit(X_train, y_train)

    y_hat = regr.predict(X_test)

    r2 = r2_score(y_test, y_hat)
    rmse = np.sqrt(np.mean((y_test - y_hat)**2))
    pearson = np.corrcoef(y_test, y_hat)[0, 1]

    fig = px.scatter(x=y_test, y=y_hat, labels={'x': 'True binding affinity', 'y': 'Predicted binding affinity'}, title=f'True vs Predicted binding affinity\nR^2 = {r2_score(y_test, y_hat):.2f}, n = {len(y_test)}, RMSE: {np.sqrt(np.mean((y_test - y_hat)**2)):.2f}, Pearson: {np.corrcoef(y_test, y_hat)[0, 1]:.2f}')
    fig.write_html(f'additional_results/GBR_strong_ablated_results_{i}.html')

    return {'r2': r2, 'rmse': rmse, 'pearson': pearson}


with Pool(cpu_count()- 4) as p:
    stats = p.map(trial, range(100))

# %%
rmses = [s['rmse'] for s in stats]
r2s = [s['r2'] for s in stats]
pearsons = [s['pearson'] for s in stats]

print(f'RMSE: {np.mean(rmses):.4f} +/- {np.std(rmses):.4f}')
print(f'R2: {np.mean(r2s):.4f} +/- {np.std(r2s):.4f}')
print(f'Pearson: {np.mean(pearsons):.4f} +/- {np.std(pearsons):.4f}')

# %% [markdown]
# ### Start from 77-dimensional vector
# 
# At each iteration, ablate the least helpful element while keeping track of which index was thrown away

# %%
from tqdm import tqdm

def train_and_test(observations, binding_affinities, num_runs=10):
    stats = []

    for i in range(num_runs):
        X_train, X_test, y_train, y_test = train_test_split(observations, binding_affinities, test_size=0.1, random_state=i)
        regr = GradientBoostingRegressor(random_state=2*i)
        regr.fit(X_train, y_train)

        y_hat = regr.predict(X_test)

        r2 = r2_score(y_test, y_hat)
        rmse = np.sqrt(np.mean((y_test - y_hat)**2))
        pearson = np.corrcoef(y_test, y_hat)[0, 1]

        stats.append({'r2': r2, 'rmse': rmse, 'pearson': pearson})

    rmses = [s['rmse'] for s in stats]
    r2s = [s['r2'] for s in stats]
    pearsons = [s['pearson'] for s in stats]

    results = {
        'rmse': {
            'mean': np.mean(rmses),
            'std': np.std(rmses)
        },
        'r2': {
            'mean': np.mean(r2s),
            'std': np.std(r2s)
        },
        'pearson': {
            'mean': np.mean(pearsons),
            'std': np.std(pearsons)
        }
    }

    return results

from multiprocessing import Pool, cpu_count

# This takes a while, but I saved the results into a pickle file to load later

print(f'Using {cpu_count() - 4} cores')

# starting with highly_selected_observations
# 10 runs for each ablated element

assert highly_selected_observations.shape[1] == len(highly_selected_observations_indices)

ablation_iter_results = []  # Keeps track of results for each iteration of ablation

# ablate one element at a time
for iter in range(highly_selected_observations.shape[1] - 1):
    print(f'Starting iteration {iter}')

    def ablate_run(index):
        ablated_observations = np.delete(highly_selected_observations, index, axis=1)
        ablation_result = train_and_test(ablated_observations, binding_affinities, num_runs=50)
        return ablation_result

    with Pool(cpu_count() - 4) as p:
        ablation_results = p.map(ablate_run, range(highly_selected_observations.shape[1]))

    with open(f'additional_results/ablation_results_{iter}.pkl', 'wb') as f:
        pickle.dump(ablation_results, f)

    # Find index with lowest RMSE (and thus the least important feature)
    min_rmse_index = np.argmin([r['rmse']['mean'] for r in ablation_results])
    min_rmse = ablation_results[min_rmse_index]['rmse']['mean']
    min_rmse_std = ablation_results[min_rmse_index]['rmse']['std']

    print(f'Finished iteration {iter}, deleted index {min_rmse_index} (original index {highly_selected_observations_indices[min_rmse_index]}), new RMSE: {min_rmse:.4f} +/- {min_rmse_std:.4f}')

    # Delete that index
    highly_selected_observations = np.delete(highly_selected_observations, min_rmse_index, axis=1)
    highly_selected_observations_indices = np.delete(highly_selected_observations_indices, min_rmse_index)

    ablation_iter_results.append({
        'min_rmse_index': min_rmse_index,
        'iteration_results': ablation_results[min_rmse_index],
        'remaining_indices': highly_selected_observations_indices
    })

    with open(f'ablation_iter_results.pkl', 'wb') as f:
        pickle.dump(ablation_iter_results, f)


# with open(f'ablation_iter_results.pkl', 'rb') as f:
#     ablation_iter_results = pickle.load(f)