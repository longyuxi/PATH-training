1. Download PDBBind refined set

```bash
cd [some_folder]
bash download_pdbbind_and_clean.sh
```

2. Modify the parameters on top of `dispatch_jobs.py` and run `python dispatch_jobs.py` to calculate the persistence diagram of each protein-ligand complex.

3. Use `make_persistence_images.ipynb` to create persistence images. *Optionally, you can use `make_persistence_images_bandwidth_1.ipynb` to create persistence images using a larger standard deviation for the Gaussians. This larger standard deviation has been found to be ineffective.*
