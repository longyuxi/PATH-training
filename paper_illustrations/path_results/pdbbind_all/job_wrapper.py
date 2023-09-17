import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent / '..'))
sys.path.append(str(Path(__file__).parent / '..' / '..'))
sys.path.append(str(Path(__file__).parent / '..' / '..' / '..'))
sys.path.append(str(Path(__file__).parent / '..' / '..' / '..' / 'gbr' / 'gbr_inference'))
from db import DB

import argparse
import traceback
from tqdm import tqdm
import inference

def job(key):
    # The part where the job actually runs, given df and idx as input

    # Each entry in the redis database should be a dictionary in the following form

    d = DB.hgetall(key)
    print(d)

    protein_file = d['protein_file']
    ligand_file = d['ligand_file']

    predictions = inference.predict(protein_file, ligand_file)
    fingerprint = inference.get_fingerprint(protein_file, ligand_file)

    # Save results
    DB.hset(key, mapping={
        **d,
        'fingerprint': str(fingerprint),
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
