import os
import pathlib
import logging
import redis
from tqdm import tqdm
import time
from ph import INDEX

logging.basicConfig(level=logging.INFO)

###############################
# Platform specific variables #
#                             #
# Change to fit to job        #
###############################


KEY_PREFIX = 'gbrph_' # Prefix of every job, as appears in the Redis database
CLUSTER = 'CS' # or 'DCC'

if CLUSTER == 'CS':
    cwd = pathlib.Path(__file__).parent.resolve()
    NUM_JOBS_TO_SUBMIT = 6000
    PYTHON_EXECUTABLE = '/usr/project/dlab/Users/jaden/mambaforge/envs/tnet2017/bin/python'
    ROOT_DIR = '/usr/project/dlab/Users/jaden/gbr-tnet/persistence'
    os.system(f'mkdir -p {ROOT_DIR}/slurm-outs')
    SBATCH_TEMPLATE = f"""#!/bin/zsh
#SBATCH --requeue
#SBATCH --chdir={ROOT_DIR}
#SBATCH --output={ROOT_DIR}/slurm-outs/%x-%j-slurm.out
#SBATCH --mem=5000M
#SBATCH --partition=compsci

source ~/.zshrc
date
hostname
conda activate tnet2017
cd {ROOT_DIR}


    """

    DB = redis.Redis(host='cybermen.cs.duke.edu', port=6379, decode_responses=True, password="topology")

    PDBBIND_BASE_FOLDER = '/usr/project/dlab/Users/jaden/pdbbind/refined-set'
    OUTPUT_BASE_FOLDER = '/usr/project/dlab/Users/jaden/gbr-tnet/persistence/persistence_diagrams'
    os.system(f'mkdir -p {OUTPUT_BASE_FOLDER}')

    ADDITIONAL_SAVE_FOLDER = ROOT_DIR + '/additional_results'
    os.system(f'mkdir -p {ADDITIONAL_SAVE_FOLDER}')

else:
    raise NotImplementedError     # Incorrect specification of cluster variable



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
        time.sleep(1)
        os.system('rm run.sh')


def populate_db(rebuild_all_keys=False):
    logging.info('Populating database')
    n_jobs = NUM_JOBS_TO_SUBMIT
    keys = [KEY_PREFIX + INDEX.at[i, 'PDB code'] for i in range(len(INDEX))]

    database_keys = DB.keys()

    for i in tqdm(range(len(INDEX))):
        k = keys[i]

        if not rebuild_all_keys and k in database_keys:
            logging.debug(f"Key {k} already exists in database")
            continue

        if rebuild_all_keys:
            DB.delete(k)

        # Get the index dataframe info as a dictionary
        info = INDEX.iloc[i].to_dict()

        db_entry = {
            'attempted': 'False',
            'error': 'False',
            'finished': 'False',
            'protein_file': f'{PDBBIND_BASE_FOLDER}/{info["PDB code"]}/{info["PDB code"]}_protein.pdb',
            'ligand_file': f'{PDBBIND_BASE_FOLDER}/{info["PDB code"]}/{info["PDB code"]}_ligand.mol2',
            'save_folder': str(OUTPUT_BASE_FOLDER),
            **info
        }


        DB.hset(k, mapping=db_entry)


    print(DB.hgetall(keys[0]))


if __name__ == '__main__':
    # rebuild_db()
    main(dry_run=True, rebuild_all_keys=True)
    main(dry_run=False)
