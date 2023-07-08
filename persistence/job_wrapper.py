import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent / '..'/ 'ph'))

import argparse
import pandas as pd
import traceback
import pickle

import dispatch_jobs
import ph
from ph import INDEX
from dispatch_jobs import DB

def job(key):
    # The part where the job actually runs, given df and idx as input

    # Each entry in the redis database should be a dictionary in the following form
    # {'attempted': 'False', 'error': 'False', 'finished': 'False', 'protein_file': '/hpc/group/donald/yl708/pdbbind/refined-set/3g2y/3g2y_protein.pdb', 'ligand_file': '/hpc/group/donald/yl708/pdbbind/refined-set/3g2y/3g2y_ligand.mol2', 'save_folder': '/hpc/group/donald/yl708/persistence-diagrams', 'PDB code': '3g2y', 'resolution': 1.31, 'release year': 2009, '-logKd/Ki': 2.0, 'Kd/Ki': 'Ki=10mM', 'reference': '3g2y.pdf', 'ligand name': 'GF4'}

    d = DB.hgetall(key)
    protein_file = d['protein_file']
    ligand_file = d['ligand_file']

    pw_opposition_diagrams = ph.get_pairwise_opposition_persistence_diagrams(protein_file, ligand_file)
    other_persistence_diagrams = ph.get_2345_persistence_diagrams(protein_file, ligand_file)

    save_file = d['save_folder'] + '/' + d['PDB code'] + '.pkl'
    d['save_file'] = save_file
    DB.hset(key, mapping=d)

    with open(save_file, 'wb') as f:
        pickle.dump({
            'pw_opposition_diagrams': pw_opposition_diagrams,
            'other_persistence_diagrams': other_persistence_diagrams
            }, f)


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
