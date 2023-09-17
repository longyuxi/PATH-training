import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent / '..'))
sys.path.append(str(Path(__file__).parent / '..' / '..'))
sys.path.append(str(Path(__file__).parent / '..' / '..' / '..'))
sys.path.append(str(Path(__file__).parent / '..' / '..' / '..' / 'persistence'))
from db import DB
from preprocessing import get_mol2_coordinates

import argparse
import traceback
from tqdm import tqdm
import inference

import analyze_gbr_importance

def job(key):
    # The part where the job actually runs, given df and idx as input

    # Each entry in the redis database should be a dictionary in the following form

    d = DB.hgetall(key)
    print(d)

    protein_file = d['Protein File']
    ligand_file = d['Ligand File']

    predictions = inference.predict(protein_file, ligand_file)
    fingerprint = inference.get_fingerprint(protein_file, ligand_file)

    num_protein_atoms = len(get_mol2_coordinates(d['Protein File']))
    num_ligand_atoms = len(get_mol2_coordinates(d['Ligand File']))
    total_atoms = num_protein_atoms + num_ligand_atoms

    # Save results
    DB.hset(key, mapping={
        **d,
        'fingerprint': str(fingerprint),
        'num_protein_atoms': str(num_protein_atoms),
        'num_ligand_atoms': str(num_ligand_atoms),
        'total_atoms': str(total_atoms),
        **predictions
    })


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
