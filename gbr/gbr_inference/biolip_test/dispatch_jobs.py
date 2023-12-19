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

import os
import pathlib
import logging
from tqdm import tqdm
import time
from pathlib import Path
import argparse
import logging

parser = argparse.ArgumentParser()
parser.add_argument( '-log',
                     '--loglevel',
                     default='warning',
                     help='Provide logging level. Example --loglevel debug, default=warning' )

args = parser.parse_args()

logging.basicConfig( level=args.loglevel.upper() )
logging.info( 'Logging now setup.' )

import sys
sys.path.append(str(Path(__file__).parent / '..' / '..' / '..'))
from db import DB

###############################
# Platform specific variables #
#                             #
# Change to fit to job        #
###############################


KEY_PREFIX = 'biolip_test_' # Prefix of every job, as appears in the Redis database
CLUSTER = 'CS' # or 'DCC'


if CLUSTER == 'CS':
    cwd = pathlib.Path(__file__).parent.resolve()
    NUM_JOBS_TO_SUBMIT = 50047
    PYTHON_EXECUTABLE = '/usr/project/dlab/Users/jaden/mambaforge/envs/tnet2017/bin/python'
    ROOT_DIR = str(Path(__file__).parent.resolve())
    os.system(f'mkdir -p {ROOT_DIR}/slurm-outs')
    SBATCH_TEMPLATE = f"""#!/bin/zsh
#SBATCH --requeue
#SBATCH --chdir={ROOT_DIR}
#SBATCH --output={ROOT_DIR}/slurm-outs/%x-%j-slurm.out
#SBATCH --mem=16000M
#SBATCH --cpus-per-task=1
#SBATCH --partition=compsci
#SBATCH --exclude=linux[1-10]

source ~/.zshrc
date
hostname
conda activate tnet2017
cd {ROOT_DIR}


    """


    ADDITIONAL_SAVE_FOLDER = ROOT_DIR + '/additional_results'
    os.system(f'mkdir -p {ADDITIONAL_SAVE_FOLDER}')

elif CLUSTER == 'DCC':
    raise NotImplementedError
    NUM_JOBS_TO_SUBMIT = 81000
    PYTHON_EXECUTABLE = '/hpc/group/donald/yl708/mambaforge/envs/tnet2017/bin/python'
    ROOT_DIR = '/hpc/group/donald/yl708/TopologyNet-2017/perturbations'
    SBATCH_TEMPLATE = f"""#!/bin/bash
#SBATCH --partition=common-old,scavenger
#SBATCH --requeue
#SBATCH --chdir={ROOT_DIR}
#SBATCH --output={ROOT_DIR}/slurm_logs/%x-%j-slurm.out
#SBATCH --mem=2500M

source ~/.bashrc
source ~/.bash_profile
date
hostname
conda activate tnet2017
cd {ROOT_DIR}


    """

    ORIGINAL_PROTEIN_FILE = '/hpc/group/donald/yl708/perturbations-data/1a4k_protein.pdb'
    LIGAND_FILE = '/hpc/group/donald/yl708/perturbations-data/1a4k_ligand.mol2'
    PERTURBATION_SAVE_FOLDER = '/hpc/group/donald/yl708/perturbations-data/'

else:
    raise Exception     # Incorrect specification of cluster variable



#############################
# Pre-execution Tests       #
#############################

# Database connection
DB.set('connection-test', '123')
if DB.get('connection-test') == '123':
    DB.delete('abc')
    logging.info('Database connection successful')
else:
    raise Exception     # Database connection failed


#############################
# Actual logic              #
#############################

# Each entry in the redis database should contain the additional information

# Run index
# Run seed
# R^2
# MSE
# Pearson
# Predictions (filename)
# Scatter plot (filename)


def main(dry_run=False, rebuild_all_keys=False):
    # Initialize database on first run
    if dry_run:
        populate_db(rebuild_all_keys=rebuild_all_keys)

    # Then submit jobs until either running out of entries or running out of number of jobs to submit
    i = 0

    database_keys = DB.keys(KEY_PREFIX + '*')
    for key in database_keys:
        if i == NUM_JOBS_TO_SUBMIT:
            break
        info = DB.hgetall(key)

        if info['finished'] == 'True' and info['error'] == 'False':
        # if info['attempted'] == 'True':
            continue
        else:
            i += 1
            # submit job for it
            if not dry_run:
                info['attempted'] = 'True'
                DB.hset(key, mapping=info)

                # sbatch run job wrapper
                sbatch_cmd = SBATCH_TEMPLATE + f'\n{PYTHON_EXECUTABLE} {str(pathlib.Path(__file__).parent) + "/job_wrapper.py"} --key {key}'

                # print(sbatch_cmd)
                with open('run.sh', 'w') as f:
                    f.write(sbatch_cmd)

                os.system(f'sbatch --job-name={key} run.sh')

    if dry_run:
        print(f'Number of jobs that would be submitted: {i}')
        time.sleep(5)
    else:
        print(f'Number of jobs submitted: {i}')


def populate_db(rebuild_all_keys=False):

    import pandas as pd

    df = pd.read_csv('./binding_affinity_data.txt', sep='\t', names=["PDB ID", "Receptor chain", "Resolution", "Binding site number code", "Ligand ID in the Chemical Component Dictionary (CCD) used by the PDB database", "Ligand chain", "Ligand serial number", "Binding site residues (with PDB residue numbering)", "Binding site residues (with residue re-numbered starting from 1)", "Catalytic site residues (different sites are separated by ';') (with PDB residue numbering)", "Catalytic site residues (different sites are separated by ';') (with residue re-numbered starting from 1)", "EC number", "GO terms", "Binding affinity by manual survey of the original literature", "Binding affinity provided by the Binding MOAD database", "Binding affinity provided by the PDBbind-CN database", "Binding affinity provided by the BindingDB database", "UniProt ID", "PubMed ID", "Residue sequence number of the ligand (field _atom_site.auth_seq_id in PDBx/mmCIF format)", "Receptor sequence"])

    from functional import seq
    import os

    BIOLIP_FOLDER = '/usr/project/dlab/Users/jaden/BioLiP_updated_set'

    # e.g. receptor/101mA.pdb
    proteins = seq(df.iterrows()).map(lambda x: 'receptor/' + x[1]['PDB ID'] + x[1]['Receptor chain'] + '.pdb').list()
    # e.g. receptor/104m_HEM_A_1.pdb
    ligands = seq(df.iterrows()).map(lambda x: 'ligand/' + x[1]['PDB ID'] + '_' + x[1]['Ligand ID in the Chemical Component Dictionary (CCD) used by the PDB database'] + '_' + x[1]['Ligand chain'] + '_' + str(x[1]['Ligand serial number']) + '.pdb').list()

    proteins = seq(proteins).map(lambda x: os.path.join(BIOLIP_FOLDER, x)).list()
    ligands = seq(ligands).map(lambda x: os.path.join(BIOLIP_FOLDER, x)).list()

    df['Protein File'] = proteins
    df['Ligand File'] = ligands

    logging.info('Populating database')
    n_jobs = NUM_JOBS_TO_SUBMIT
    keys = [KEY_PREFIX + str(i) for i in range(n_jobs)]

    database_keys = DB.keys()

    for k in tqdm(keys):
        if not rebuild_all_keys and k in database_keys:
            logging.debug(f"Key {k} already exists in database")
            continue

        if rebuild_all_keys:
            DB.delete(k)

        df_index = int(k[len(KEY_PREFIX):])

        # Add the entry to the database
        DB.hset(k, mapping={
            'attempted': 'False',
            'error': 'False',
            'finished': 'False',
            **df.iloc[df_index].to_dict()
        })

    print(f'Example entry: {DB.hgetall(keys[0])}')


def rebuild_db():
    raise NotImplementedError

def get_db():
    # Pinnacle of OOP
    return DB

if __name__ == '__main__':
    # rebuild_db()
    main(dry_run=True, rebuild_all_keys=True)
    main(dry_run=False)
