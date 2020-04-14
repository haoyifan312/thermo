import numpy as np


def normalize_array(array):
    """
    normalize a vector with sum to 1
    """
    size = array.size
    if size == 0:
        raise Exception("normalize array of size 0")
    sum = np.sum(array)
    if sum == 0.0:
        raise Exception("normalize array with sum 0")
    return array/sum
