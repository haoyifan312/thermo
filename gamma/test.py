# test UNIQUAC, UNIFAC code
from unifac import *
from uniquac import *
import numpy as np

# water
# ethanol
# hexane

inflow = ['water', 'ethanol', 'hexane']
water_group = np.array([[16, 1]])
ethanol_group = np.array([[1, 1], [2, 1], [14, 1]])
hexane_group = np.array([[1, 2], [2, 4]])
inflow_group = [water_group, ethanol_group, hexane_group]
rq = get_r_q_list(inflow_group)
size = len(inflow_group)
n = np.random.rand(size)
x = n/np.sum(n)

gamC = staverman_guggenheim(rq, x)

group_set = get_group_set(inflow_group)

end =1
