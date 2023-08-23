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
