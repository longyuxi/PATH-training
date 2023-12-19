## This file is part of PATH, which is part of OSPREY 3.0
##
## OSPREY Protein Redesign Software Version 3.0
## Copyright (C) 2001-2023 Bruce Donald Lab, Duke University
##
## OSPREY is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License version 2
## as published by the Free Software Foundation.
##
## You should have received a copy of the GNU General Public License
## along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
##
## OSPREY relies on grants for its development, and since visibility
## in the scientific literature is essential for our success, we
## ask that users of OSPREY cite our papers. See the CITING_OSPREY
## document in this distribution for more information.
##
## Contact Info:
##    Bruce Donald
##    Duke University
##    Department of Computer Science
##    Levine Science Research Center (LSRC)
##    Durham
##    NC 27708-0129
##    USA
##    e-mail: www.cs.duke.edu/brd/
##
## <signature of Bruce Donald>, Mar 1, 2023
## Bruce Donald, Professor of Computer Science

import sys
from pathlib import Path
sys.path.append(str((Path(__file__)/ '..' / '..' / '..' / 'persistence').resolve()))
import pickle

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import r2_score
import plotly.express as px
from joblib import Memory
memory = Memory("__joblib_cache__", verbose=0)

import ph
from gudhi.representations.vector_methods import PersistenceImage

highly_selected_observations_indices = [16829, 16930, 17031, 17818, 17819, 17917, 17918, 17919, 17920, 18019, 18020, 18021, 18115, 18116, 18117, 18120, 18121, 18214, 18215, 18216, 18217, 18218, 18313, 18314, 18315, 18316, 18317, 18318, 18413, 18415, 18416, 18417, 18418, 18514, 18516, 18517, 18610, 18615, 18710, 18711, 18811, 18812, 18908, 18912, 18913, 19110, 287524, 287620, 287719, 287720, 287721, 287818, 287821, 287822, 287919, 287921, 287922, 287923, 288020, 288022, 288121, 288219, 288314, 288415, 288516, 31500, 31700, 31701, 31702, 61301, 91000, 91001, 91002, 91802, 92100, 92101, 92102]

tiny_feature_subset_indices = [17819,  18318,  18514,  18908, 288020, 288516,  31500,  31702, 61301,  91001]

# Loads regressors if they are present, otherwise trains them on PDBBind dataset and saves them
def train_or_load_regrs(force_retrain=False, regr_77_path='/usr/project/dlab/Users/jaden/gbr-tnet/gbr/gbr_inference/regr_with_77_elements.pkl', regr_10_path='/usr/project/dlab/Users/jaden/gbr-tnet/gbr/gbr_inference/regr_with_10_elements.pkl', optimal_mini_gbr_path='/usr/project/dlab/Users/jaden/gbr-tnet/gbr/gbr_inference/optimal_mini_gbr.pkl'):
    all_gbrs_exist = Path(regr_77_path).exists() and Path(regr_10_path).exists() and Path(optimal_mini_gbr_path).exists()

    if all_gbrs_exist and not force_retrain:
        with open(regr_77_path, 'rb') as f:
            regr_with_77_elements = pickle.load(f)
        with open(regr_10_path, 'rb') as f:
            regr_with_10_elements = pickle.load(f)
        with open(optimal_mini_gbr_path, 'rb') as f:
            optimal_mini_gbr = pickle.load(f)

        print(f'Loaded from {regr_77_path}, {regr_10_path} and {optimal_mini_gbr_path}')

        return regr_with_77_elements, regr_with_10_elements, optimal_mini_gbr

    else:
        print('Training GBR')
        num_datapoints_per_file = 500

        # Load the training data
        observations = []
        binding_affinities = []

        NPYS_FOLDER = Path('/usr/project/dlab/Users/jaden/gbr-tnet/persistence/persistence_images')

        for file_index in range(len(list(NPYS_FOLDER.glob('binding_affinities_*.npy')))):
            with open(NPYS_FOLDER / f'observations_{file_index}.npy', 'rb') as f:
                observations.append(np.load(f))

            with open(NPYS_FOLDER / f'binding_affinities_{file_index}.npy', 'rb') as f:
                binding_affinities.append(np.load(f))

        observations = np.concatenate(observations)
        binding_affinities = np.concatenate(binding_affinities)

        print(observations.shape)
        print(binding_affinities.shape)

        observations = observations.reshape(observations.shape[0], -1)
        observation_indices = np.arange(observations.shape[1])

        print(observations.shape)
        print(binding_affinities.shape)

        highly_selected_observations = observations[:, highly_selected_observations_indices]

        # Train on the entire dataset
        regr_with_77_elements = GradientBoostingRegressor()
        regr_with_77_elements.fit(highly_selected_observations, binding_affinities)


        tiny_feature_subset = observations[:, tiny_feature_subset_indices]
        regr_with_10_elements = GradientBoostingRegressor()
        regr_with_10_elements.fit(tiny_feature_subset, binding_affinities)

        optimal_gbr_kwargs = {'max_depth': 3, 'learning_rate': 0.4, 'n_estimators': 13}
        optimal_mini_gbr = GradientBoostingRegressor(**optimal_gbr_kwargs)
        optimal_mini_gbr.fit(tiny_feature_subset, binding_affinities)

        with open(regr_77_path, 'wb') as f:
            pickle.dump(regr_with_77_elements, f)

        with open(regr_10_path, 'wb') as f:
            pickle.dump(regr_with_10_elements, f)
        
        with open(optimal_mini_gbr_path, 'wb') as f:
            pickle.dump(optimal_mini_gbr, f)

        return regr_with_77_elements, regr_with_10_elements, optimal_mini_gbr

regr_with_77_elements, regr_with_10_elements, optimal_mini_gbr = train_or_load_regrs()


def diagram_to_image(diagram):
    if diagram is None:
        return np.zeros((3, 10000))
    else:
        pim = PersistenceImage(bandwidth=0.2355, resolution=[100,100], im_range=[0, 50, 0, 50])

        # check that the third column is indeed only 0, 1, 2
        # assert (lambda unique_elements_in_array=np.unique(diagram[0, :, 2])
                    # : np.array_equal(unique_elements_in_array, np.array([0, 1, 2])))()

        # Now reshape into (3, n, 2)
        def get_nd_diagram(diagram, dim):
            filtered_diagrams = list(filter(lambda x: x[2] == dim, diagram[0, :, :]))
            return np.array(filtered_diagrams)[:, :2]

        gudhi_format_diagrams = [get_nd_diagram(diagram, 0), get_nd_diagram(diagram, 1), get_nd_diagram(diagram, 2)]


        diagrams_clipped = [np.clip(diagram, 0, 50) for diagram in gudhi_format_diagrams]
        imgs = pim.fit_transform(diagrams_clipped)

        return imgs


# Takes protein file, ligand file, and predicts
@memory.cache
def get_all_images(protein_file, ligand_file):
    # do persistent homology on protein_file and ligand_file
    pw_opposition_diagrams = ph.get_pairwise_opposition_persistence_diagrams(protein_file, ligand_file)
    other_persistence_diagrams = ph.get_2345_persistence_diagrams(protein_file, ligand_file)

    all_diagrams = pw_opposition_diagrams + other_persistence_diagrams
    all_images = list(map(diagram_to_image, all_diagrams))
    return all_images

def predict(protein_file, ligand_file):
    print('Computing persistence images')
    all_images = get_all_images(protein_file, ligand_file)

    print('Predicting')
    observations = np.array(all_images).flatten().reshape(1, -1)
    observation_77 = observations[:, highly_selected_observations_indices]
    observation_10 = observations[:, tiny_feature_subset_indices]

    prediction_77 = regr_with_77_elements.predict(observation_77)[0]
    prediction_10 = regr_with_10_elements.predict(observation_10)[0]
    prediction_mini = optimal_mini_gbr.predict(observation_10)[0]

    return {
        'prediction_77': prediction_77,
        'prediction_10': prediction_10,
        'prediction_mini': prediction_mini
    }

def get_fingerprint(protein_file, ligand_file):
    all_images = get_all_images(protein_file, ligand_file)
    observations = np.array(all_images).flatten().reshape(1, -1)
    observation_10 = observations[:, tiny_feature_subset_indices]
    return observation_10[0]


if __name__ == '__main__':
    protein_file = '/usr/project/dlab/Users/jaden/gbr-tnet/example_pdbbind_folder/1a1e/1a1e_protein.pdb'
    ligand_file = '/usr/project/dlab/Users/jaden/gbr-tnet/example_pdbbind_folder/1a1e/1a1e_ligand.mol2'
    res = predict(protein_file, ligand_file)
    print(res)