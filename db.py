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

## Constants for database and index location

import redis
from pathlib import Path
from persistence.preprocessing import load_pdbbind_data_index

CLUSTER = 'CS' # or 'DCC'


if CLUSTER == 'CS':
    DB = redis.Redis(host='cybermen.cs.duke.edu', port=6379, decode_responses=True, password="topology")
elif CLUSTER == 'DCC':
    DB = redis.Redis(host='dcc-login-03', port=6379, decode_responses=True, password="topology")
else:
    raise Exception     # Incorrect specification of cluster variable


if CLUSTER == 'local':
    INDEX_LOCATION = Path('/home/longyuxi/Documents/mount/pdbbind-dataset/index/INDEX_refined_data.2020')
    BASE_FOLDER = Path('/home/longyuxi/Documents/mount/pdbbind-dataset/refined-set')
elif CLUSTER == 'DCC':
    INDEX_LOCATION = Path('/hpc/group/donald/yl708/pdbbind/index/INDEX_refined_data.2020')
    BASE_FOLDER = Path('/hpc/group/donald/yl708/pdbbind/refined-set')
elif CLUSTER == 'CS':
    INDEX_LOCATION = Path('/usr/project/dlab/Users/jaden/pdbbind/index/INDEX_refined_data.2020')
    BASE_FOLDER = Path('/usr/project/dlab/Users/jaden/pdbbind/refined-set')
else:
    raise NotImplementedError


INDEX = load_pdbbind_data_index(INDEX_LOCATION)

PDBBIND_BASE_FOLDER = '/usr/project/dlab/Users/jaden/pdbbind/refined-set'
