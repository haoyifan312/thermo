# import sys, os
# cwd = os.getcwd()
# sys.path.append(os.path.join(cwd,'gamma'))
# sys.path.append(os.path.join(cwd,'test'))

# test UNIQUAC, UNIFAC code
from unifac import get_group_set, rq_from_db
from uniquac import staverman_guggenheim
from testHelper import *
import numpy as np

# water
# ethanol
# hexane

inflow = ['water', 'ethanol', 'hexane']
water_group = np.array([[16, 1]])
ethanol_group = np.array([[1, 1], [2, 1], [14, 1]])
hexane_group = np.array([[1, 2], [2, 4]])
inflow_group = [water_group, ethanol_group, hexane_group]
rq = rq_from_db(inflow_group)


def test_rq():
    rq_save = np.array([[0.92, 2.5755, 4.4998], [1.4, 2.588, 3.856]])
    assert verify_vec_rel(rq, rq_save)


size = len(inflow_group)
n = np.ones(size)
x = n/np.sum(n)
gamC = staverman_guggenheim(rq[0, :], rq[1, :], x)


def test_gamC():
    gamC_save = np.array([0.17683701, 0.00310259, 0.02007396])
    assert verify_vec_rel(gamC, gamC_save)


group_set = get_group_set(inflow_group)


def test_groupset():
    group_set_save = np.array([16,  1,  2, 14])
    assert verify_vec_rel(group_set, group_set_save)


end = 1
