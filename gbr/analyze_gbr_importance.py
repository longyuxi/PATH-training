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

import pickle
from sklearn.inspection import permutation_importance
import numpy as np
import plotly.express as px
import pandas as pd
import numpy as np
import pickle


# Permutation importance takes an uncanny amount of time to run for the 1.2M persistent homology features, and the impurity-based feature importance seems already sufficient, so we will not use permutation importance.
def analyze(regr, seed, impurity_importance_html, sorted_impurity_importances_pckl, permutation_importance_html, sorted_permutation_importances_pckl, test_observations, test_binding_affinities):
    analyze_impurity_importance(regr, impurity_importance_html, sorted_impurity_importances_pckl)
    # analyze_permutation_importance(regr, seed, permutation_importance_html, sorted_permutation_importances_pckl, test_observations, test_binding_affinities)


def analyze_impurity_importance(regr, impurity_importance_html, sorted_impurity_importances_pckl):
    feature_importance = regr.feature_importances_

    top_feature_threshold = 0.005
    num_top_features = len(list(filter(lambda x: x > top_feature_threshold, feature_importance)))

    df = pd.DataFrame({'feature': np.arange(len(feature_importance)), 'importance': feature_importance})

    # Sort by importance
    df = df.sort_values('importance', ascending=False)

    # Get top features
    df_top = df[:num_top_features]

    # Convert feature column to string
    df_top['feature'] = df_top['feature'].astype(str)

    # Plotly bar plot
    fig = px.bar(df_top, x='feature', y='importance', title=f'Features with importance > {top_feature_threshold} found via <b>impurity based feature importance</b>')

    fig.write_html(impurity_importance_html)

    with open(sorted_impurity_importances_pckl, 'wb') as f:
        pickle.dump(df, f)

    # Sort by importance

    df = df.sort_values('importance', ascending=True)

    # Add cumulative importance column
    df['cumulative_importance'] = df['importance'].cumsum()

    # Area under Lorenz curve
    auc = df['importance'].mean()
    print(f'Area under Lorenz curve: {auc}')


def analyze_permutation_importance(regr, seed, permutation_importance_html, sorted_permutation_importances_pckl, test_observations, test_binding_affinities):
    # ### Permutation importance
    # Which is a better method according to https://explained.ai/rf-importance/
    # Load test set

    from sklearn.model_selection import train_test_split

    test_observations = test_observations.reshape(test_observations.shape[0], -1)  # Flatten out for SVM

    X_test, y_test = test_observations, test_binding_affinities

    from sklearn.inspection import permutation_importance

    permutation_result = permutation_importance(
        regr, X_test, y_test, n_repeats=10, random_state=seed, n_jobs=-1,
    )

    permutation_top_feature_threshold = 0.005
    num_top_permutation_features = len(list(filter(lambda x: x > permutation_top_feature_threshold, permutation_result.importances_mean)))

    df = pd.DataFrame({'feature': np.arange(len(permutation_result.importances_mean)), 'permutation_importance': permutation_result.importances_mean})

    # Sort by importance
    df = df.sort_values('permutation_importance', ascending=False)

    # Get top features
    df_top = df[:num_top_permutation_features]

    # Convert feature column to string
    df_top['feature'] = df_top['feature'].astype(str)

    # Plotly bar plot
    fig = px.bar(df_top, x='feature', y='permutation_importance', title=f'Features with importance > {permutation_top_feature_threshold} found via <b>permutation importance</b>')

    fig.write_html(permutation_importance_html)

    with open(sorted_permutation_importances_pckl, 'wb') as f:
        pickle.dump(df, f)

