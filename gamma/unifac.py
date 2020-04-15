import numpy as np
import pandas as pd


def get_r_q_list(inflow_group):
    """
    Get list of UNIQUAC r and q parameters for the list of inflow species
    """
    groupdata_csv = r'.\data\unifac\group_data.csv'
    groupdata = pd.read_csv(groupdata_csv, index_col='subgroup')
    spsize = len(inflow_group)
    rq = np.zeros((2, spsize))
    for i, molecule in enumerate(inflow_group):
        for group in molecule:
            subgroup = group[0]
            ngroup = group[1]
            rq[0, i] += groupdata.at[subgroup, 'R']*ngroup
            rq[1, i] += groupdata.at[subgroup, 'Q']*ngroup
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
