# Predicting Affinity Through Homology (PATH): Interpretable Binding Affinity Prediction with Persistent Homology

This is the full repository containing scripts for feature construction, feature selection, training, and testing of the interpretable protein-ligand binding affinity prediction algorithm *Predicting Affinity Through Homology (PATH)*.

The inference code for PATH can be found [in the OSPREY3 package](https://github.com/donaldlab/OSPREY3/tree/main/src/main/python/path). My implementation of TNet-BP from TopologyNet (Cang and Wei, 2017) can be found [here](https://github.com/longyuxi/TopologyNet-2017).

# Prerequisites

## 1. Python environment
1. Install [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/download.html) or [mamba](https://mamba.readthedocs.io/en/latest/mamba-installation.html) (you only need one of the two).
2. Run `conda env create -f tnet2017.yml` or `mamba env create -f tnet2017.yml` to create the `tnet2017` environment. Activate this environment by `conda activate tnet2017`.
3. `pip install -r requirements.txt` to install additional dependencies for this project that are not available through conda.

## 2. Redis server

*I execute my computational tasks (jobs) in a compute cluster managed by SLURM. To keep track of the statuses of jobs and their results, I use a Redis database and a custom MapReduce-style system.*

I first explain my custom MapReduce-style system. This system consists of two scripts, `job_wrapper.py` and `dispatch_jobs.py`, a SLURM scheduler, and a Redis database. If you are running these scripts in a SLURM cluster, you will need to modify the headers of the temporary shell scripts (see below) to fit the configuration of your cluster. If you are executing these scripts on a compute cluster with a different job scheduler, more changes will need to be made according how compute jobs are submitted on your cluster.

1. Each task is associated with a set of sequentially numbered key starting from a prefix, which is reflected in the `KEY_PREFIX` variable in `dispatch_jobs.py`.
2. `dispatch_jobs.py` will create an entry in the database for each key containing information about the job and the fields {`started`, `finished`, `error`} all set to `False`. It then submits by creating temporary shell scripts that execute `python job_wrapper.py --key {k}` and submit these shell scripts to the SLURM scheduler.
3. `job_wrapper.py` contains the instructions for execution when the work is allocated to a scheduler.

As mentioned, a Redis database is used for managing jobs submitted to the SLURM batch cluster. To set up this database,

1. Build and install Redis via [https://redis.io/docs/getting-started/installation/install-redis-from-source/].
2. Optionally, add the `src` folder of Redis to path.
3. Create a `redis.conf` file somewhere and set a default password by putting e.g. `requirepass topology` in that file.
4. Start the redis server on a host with your `redis.conf` and adjust the `DB` constant in `dispatch_jobs.py` accordingly.

## 3. Clone this repository

1. Install [git lfs](https://git-lfs.com/).
2. Clone this repository.

# Persistence image construction

The `persistence` folder contains scripts that construct persistence images using the opposition distance. Start with `persistence/README.md`

# Feature selection, model training, and validation

The `gbr` folder contains of scripts that perform feature selection, model training, and validation of the gradient boosting regressor on persistence images to eventually lead to the formulation of persistence fingerprint and PATH. Start with `gbr/README.md`.



# Citation requirements

See [CITING_OSPREY.txt](CITING_OSPREY.txt).
