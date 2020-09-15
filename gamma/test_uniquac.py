import numpy as np

import os
print(os.getcwd())

from uniquac import UNIQUAC, get_tauij
from testHelper import *

def test_uniquac():
    # hexane heptane pentane

    tau12 = [0,	0,	-14.4923,	16.0694, 0, 0, 0, 0]
    tau13 = [0,	0,	91.9378, -101.2586,	0,	0,	0,	0]
    tau23 = [-1.54345,	0.957882,	231.863,	-118.357,	0,	0,	0,	0]

    r = np.array([4.49967,	5.17403,	3.82531])
    q = np.array([3.856,	4.396,	3.316])
    tauij = np.ones((3, 3))
    tauij[0, 1] = get_tauij(350, np.array([1.0, 0.0, 0.0, 0.0, 0.0]))
    tauij[1, 0] = get_tauij(350, np.array([2.0, 0.0, 0.0, 0.0, 0.0]))
    tauij[0, 2] = get_tauij(350, np.array([0.0, 1000.0, 0.0, 0.0, 0.0]))
    tauij[2, 0] = get_tauij(350, np.array([0.0, 2000.0, 0.0, 0.0, 0.0]))
    tauij[1, 2] = get_tauij(350, np.array([0.0, 0.0, 0.1, 0.0, 0.0]))
    tauij[2, 1] = get_tauij(350, np.array([0.0, 0.0, 0.2, 0.0, 0.0]))
    uniquac = UNIQUAC(3, r, q, tauij)
    uniquac.x = np.array([.1, .2, .7])
    uniquac_save = np.array([1.03253275e-14, 1.31722851e-01, 7.76709186e-02])
    uniquac_out = np.exp(uniquac.compute())
    assert verify_vec_rel(uniquac_out, uniquac_save)

    #print(np.exp(uniquac.compute()))