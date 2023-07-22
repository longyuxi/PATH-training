import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent / '..' / 'gbr_iter_0'))

import os
import pickle
import hashlib
from functools import wraps
import inspect

import dispatch_jobs
DB = dispatch_jobs.get_db()

CACHING_DIRECTORY = Path('/usr/project/dlab/Users/jaden/_cache')

def redis_cache(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        # Generate a unique key for each function and its arguments
        print(inspect.getsource(func))
        key = 'cache_' + str(func.__name__) + str(args) + str(kwargs) + str(hashlib.md5(inspect.getsource(func).encode()).hexdigest())

        # Check if we already have a result in cache
        result_pickle_filename = DB.get(key)
        
        if result_pickle_filename is not None:
            # Load result from the pickle file
            try:
                with open(result_pickle_filename, 'rb') as f:
                    print(f'Loaded result from {result_pickle_filename}')
                    return pickle.load(f)
            except FileNotFoundError:
                # If the file doesn't exist, we need to run the function
                pass

        # If the result isn't cached, we need to run the function
        result = func(*args, **kwargs)
        
        # Pickle the result
        result_pickle_filename = str((CACHING_DIRECTORY / f'{key}.pickle').resolve())
        with open(result_pickle_filename, 'wb') as f:
            print(f'Cached results to {result_pickle_filename}')
            pickle.dump(result, f)
            
        # Store the filename in Redis
        DB.set(key, result_pickle_filename)
        
        return result
    return wrapper


@redis_cache
def dummy_function(a, b):
    import time
    print('Running dummy function')
    time.sleep(5)
    return a + b

if __name__ == '__main__':
    print(dummy_function(1, 2))
    print(dummy_function(1, 2))