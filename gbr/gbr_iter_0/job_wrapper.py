import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent / '..'))
sys.path.append(str(Path(__file__).parent / '..' / '..'))
from db import DB

import argparse
import traceback
from gtda.diagrams import PersistenceImage
import numpy as np
from tqdm import tqdm

import analyze_gbr_importance

pim = PersistenceImage()

def job(key):
    # The part where the job actually runs, given df and idx as input

    # Each entry in the redis database should be a dictionary in the following form
    # {'error': 'False', 'finished': 'False', 'pearson': '0', 'rmse': '0', 'seed': '-3981833901435230303', 'prediction_results_file': '/usr/project/dlab/Users/jaden/tnet2017-new/svm/gbr_reruns/additional_results/gbr_rerun_90_10_0_predictions.npy', 'r2': '0', 'train_ratio': '0.9', 'scatter_plot_file': '/usr/project/dlab/Users/jaden/tnet2017-new/svm/gbr_reruns/additional_results/gbr_rerun_90_10_0_scatter_plot.html', 'attempted': 'True'}

    # Need to set r2, rmse, pearson


    d = DB.hgetall(key)
    print(d)

    seed = int(d['seed']) % (2**32)
    train_ratio = float(d['train_ratio'])



    import numpy as np
    from pathlib import Path
    import pickle

    # Loading the files back in
    observations = []
    binding_affinities = []

    NPYS_FOLDER = Path('/usr/project/dlab/Users/jaden/gbr-tnet/persistence/persistence_images')
    # NPYS_FOLDER = Path('/tmp/yl708/persistence_images')

    for file_index in tqdm(range(len(list(NPYS_FOLDER.glob('binding_affinities_*.npy'))))):
        with open(NPYS_FOLDER / f'observations_{file_index}.npy', 'rb') as f:
            observations.append(np.load(f))

        with open(NPYS_FOLDER / f'binding_affinities_{file_index}.npy', 'rb') as f:
            binding_affinities.append(np.load(f))

    observations = np.concatenate(observations)
    binding_affinities = np.concatenate(binding_affinities)

    if 'test' in d:
        print('Test key detected, using only 20 samples')
        observations = observations[:20]
        binding_affinities = binding_affinities[:20]

    print(observations.shape)
    print(binding_affinities.shape)


    from sklearn.model_selection import train_test_split

    observations = observations.reshape(observations.shape[0], -1)  # Flatten out for SVM

    X_train, X_test, y_train, y_test = train_test_split(observations, binding_affinities, train_size=train_ratio, random_state=seed)

    del observations

    from sklearn.ensemble import GradientBoostingRegressor

    regr = GradientBoostingRegressor()
    print('Starting fit')

    regr.fit(X_train, y_train)

    import plotly.express as px
    from sklearn.metrics import r2_score

    print('Starting predict')

    y_hat = regr.predict(X_test)

    fig = px.scatter(x=y_test, y=y_hat, labels={'x': 'True binding affinity', 'y': 'Predicted binding affinity'}, title=f'True vs Predicted binding affinity\nR^2 = {r2_score(y_test, y_hat):.2f}, n = {len(y_test)}, RMSE: {np.sqrt(np.mean((y_test - y_hat)**2)):.2f}, Pearson: {np.corrcoef(y_test, y_hat)[0, 1]:.2f}')
    fig.write_html(d['scatter_plot_file'])

    np.save(d['prediction_results_file'], np.stack((y_test, y_hat)))

    print(f'RMSE: {np.sqrt(np.mean((y_test - y_hat) ** 2)): .4f}')
    print(f'Pearson correlation: {np.corrcoef(y_test, y_hat)[0, 1]: .4f}')

    d['r2'] = f'{r2_score(y_test, y_hat):.4f}'
    d['rmse'] = f'{np.sqrt(np.mean((y_test - y_hat) ** 2)):.4f}'
    d['pearson'] = f'{np.corrcoef(y_test, y_hat)[0, 1]:.4f}'
    DB.hset(key, mapping=d)

    analyze_gbr_importance.analyze(regr, seed=seed, impurity_importance_html=d['impurity_importance_html'], sorted_impurity_importances_pckl=d['sorted_impurity_importances_pckl'], permutation_importance_html=d['permutation_importance_html'], sorted_permutation_importances_pckl=d['sorted_permutation_importances_pckl'], test_observations=X_test, test_binding_affinities=y_test)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--key')

    args = parser.parse_args()
    key = args.key

    print('key', key)
    print('job started')
    d = DB.hgetall(key)
    d['attempted'] = 'True'
    DB.hset(key, mapping=d)

    try:
        job(key)
        print('job finished')
        d = DB.hgetall(key)
        d['finished'] = 'True'
        d['error'] = 'False'
        DB.hset(key, mapping=d)
        print('job success')

    except Exception as err:
        print(Exception, err)
        print(traceback.format_exc())
        print('job error')

        d = DB.hgetall(key)
        d['finished'] = 'True'
        d['error'] = 'True'
        DB.hset(key, mapping=d)
