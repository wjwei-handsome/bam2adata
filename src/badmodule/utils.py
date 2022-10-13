import sys
import time

from badmodule.logger import logger


def _warn(msg):
    logger.warning(f"{msg}")

def _error(msg):
    logger.error(f"{msg}")
    sys.exit(1)

def _info(msg):
    logger.info(f"{msg}")

def get_time(func):
    """ a decorator to get the time of a function and return it """
    def wrapper(*args, **kwargs):
        start = time.time()
        res = func(*args, **kwargs)
        end = time.time()
        _info(f"Process {args[0]['gene_id']}  {end - start}s ✅")
        return res
    return wrapper