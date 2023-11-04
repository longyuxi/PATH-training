# Feature and model selection

## 1.2M features → ~77 features

The scripts in `gbr_iter_0` makes 100 random trainings of the gradient boosting regressor on persistence images and note down the features with their mean decrease in impurity.

```
cd gbr_iter_0
python dispatch_jobs.py
```

*We have also tried to measure feature importance by [permutation importance](https://scikit-learn.org/stable/modules/permutation_importance.html). But the execution of permutation importance measurement took excessively long due to the size of persistence images. Therefore, this feature is disabled.*

Then `gbr_iter_1/collect_gbr_iter_0_results_glitchy.ipynb` collects results from the above run and identifies features in persistence images that were present at least once among top 10 mean decrease in impurity across the 100 runs. In our trials, we identified 77 features that satisfy this criterion.

Copy the indices of important features at the bottom of `gbr_iter_1/collect_gbr_iter_0_results_glitchy.ipynb` and paste them into the variable `highly_selected_observations_indices` in `gbr_iter_1/gbr_ablation.ipynb`.

## ~77 features → 10 features i.e., persistence fingerprint

Now we select the most important features among these ~77 features with the *ablate-and-test* procedure. Use the cells in `gbr_iter_1/gbr_ablation.ipynb` before the cell that contains `tiny_feature_subset_indices`. Read the graph that shows *RMSE vs number of remaining features* and select the smallest number of features that still preserve a good accuracy in your opinion. We found that 10 features still lead to a decent accuracy. Identify the indices of these features using the cell right below the graph. Copy this list of indices into the variable `tiny_feature_subset_indices`.

## Persistence fingerprint → PATH

The remainder of `gbr_iter_1/gbr_ablation.ipynb` identifies the optimal hyperparameters (number of trees, tree depth, and learning rate). The optimal algorithm is named PATH.

# Model testing

## Alternative regressors

Testing of alternative regressors on persistence fingerprint is done in `gbr_iter_1/regressors_test/regressors_test.ipynb`.

## Generalizability of PATH

`gbr_inference/biolip_test` contains scripts to test the generalizability of PATH on the BioLiP dataset. Instructions for downloading the BioLiP dataset are as follows:

1. Run `https://zhanggroup.org/BioLiP/download/download_all_sets.pl` to where you want to download the BioLIP structure dataset.
2. Configure the appropriate value for the variable `BIOLIP_FOLDER`.

Then in the folder `gbr_inference/biolip_test`, first modify the parameters of `dispatch_jobs.py` and run it to perform inference on the BioLiP dataset. Then use `collect_biolip_results.ipynb` to collect the inference results.

`gbr_inference` also contains scripts to test the inference capability of PATH on protein-ligand complexes measured in the [OSPREY3 paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6391056/).

