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
sys.path.append(str(Path(__file__).parent / '..'/ 'ph'))
sys.path.append(str(Path(__file__).parent / '..'))

import argparse
import traceback
import pickle

import ph
from db import DB

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
