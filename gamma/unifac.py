import numpy as np
import pandas as pd
import os, sys
import csv, sqlite3
cwd = os.getcwd()
sys.path.append(os.path.join(cwd,'util'))
from util import query_db



def rq_from_db(inflow_group):
    """
    Get list of UNIQUAC r and q parameters for the list of inflow species
    """
    # get a set of subgroups
    group_set = get_group_set(inflow_group)
    s_set = [str(i) for i in group_set]
    line = ','.join(s_set)

    #connect to db file
    cwd = os.getcwd()
    db_path = os.path.join(cwd,r'data\UNIFAC.db')
    rq_data = query_db(db_path, 'SUBGROUP, R, Q', 'GROUP_INFO', f'SUBGROUP IN ({line})')

    # convert from list to dic
    r_dic = {}
    q_dic = {}
    for subgroup, r, q in rq_data:
        r_dic[subgroup] = r
        q_dic[subgroup] = q

    spsize = len(inflow_group)
    rq = np.zeros((2, spsize))
    for i, molecule in enumerate(inflow_group):
        for group in molecule:
            subgroup = group[0]
            ngroup = group[1]
            rq[0, i] += r_dic[subgroup]*ngroup
            rq[1, i] += q_dic[subgroup]*ngroup
    return rq


def get_group_set(inflow_group):
    """
    Get the set of UNIFAC groups based on inflow species
    """
    group_list = np.empty((0), dtype=int)
    for molecule in inflow_group:
        groups = molecule[:, 0]
        group_list = np.concatenate((group_list, groups))

    return np.array(list(set(group_list)))
