import numpy as np
import pandas as pd
import os, sys
import csv, sqlite3
cwd = os.getcwd()
sys.path.append(os.path.join(cwd,'util'))
sys.path.append(os.path.join(cwd,'gamma'))
from util import query_db_one, query_db
from mdb_query import query
from thermo_common import normalize_array


def query_rq(line):
    # connect to db file
    cwd = os.getcwd()
    db_path = os.path.join(cwd, r'data\UNIFAC.db')
    # rq_data = query_db(db_path, 'SUBGROUP, R, Q', 'GROUP_INFO', f'SUBGROUP IN ({line})')
    rq_data = []
    for group in line.split(','):
        rq_data.append(query_db_one(db_path, 'SUBGROUP, R, Q', 'GROUP_INFO', f'SUBGROUP IN ({group})')[0])

    # convert from list to dic
    r_dic = {}
    q_dic = {}
    for subgroup, r, q in rq_data:
        r_dic[subgroup] = r
        q_dic[subgroup] = q
    
    return {"r_dic": r_dic, "q_dic": q_dic}

def rq_from_db(inflow_group):
    """
    Get list of UNIQUAC r and q parameters for the list of inflow species
    inflow_group - [
        [[group1, num], [group2, num], [group3, num]] - species1,
        [[group1, num], [group2, num]] - species2,
        [[group1, num], [group2, num], [group3, num]] - species3
    ]
    """
    # get a set of subgroups
    group_set = get_group_set(inflow_group)
    s_set = [str(i) for i in group_set]
    line = ','.join(s_set)

    # connect to db file
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
    group_list = np.array([group[0] for molecule in inflow_group for group in molecule])
    return np.array(list(set(group_list)))


class UNIFAC:
    def __init__(self, cas_list):
        self.cas_list = cas_list

    def update_cas_at(self, i, cas):
        """
        update cas list at index i
        """
        if i < len(self.cas_list):
            self.cas_list[i] = cas
        else:
            raise Exception('species index greater than list lenght')

    def append_cas_list(self, cas):
        """
        append cas to the species list
        """
        self.cas_list.append(cas)

    def setup_from_db(self):
        """
        setup list of group list by species from mongodb
        """
        self.isSetup = False

        # check species list size
        if len(self.cas_list) == 0:
            raise Exception('UNIFAC species list is empty')

        self.groups_by_species = query(self.cas_list) 

        self.isSetup = True

    def setup_all_species(self, vec):
        """
        setup list of group list by species from outside
        """
        self.groups_by_species = vec
        self.isSetup = True

    def update_groups_for_one(self, i, groups):
        """
        update groups of one species
        """
        if i < len(self.groups_by_species):
            self.groups_by_species[i] = groups
        else:
            raise Exception('index is greater than species list')

    def generate_group_lists(self):
        """
        # groups and counts for each species
        """
        self.subgroup_list = []
        self.maingroup_list = []
        try:
            for molecule in self.groups_by_species:
                self.subgroup_list.append([(g['subGroup'], g['count']) for g in molecule['groups']])
                self.maingroup_list.append([(g['mainGroup'], g['count']) for g in molecule['groups']])
        except Exception:
            raise Exception('species group list is empty')

    def setup_r_q(self):
        """
        self.subgroup_list/self.maingroup_list - group&count for species
        self.subgroup_set/self.maingroup_set - group set for whole system
        """

        # self.subgroup_list, self.maingroup_list
        self.generate_group_lists()

        # self.subgroup_set, self.maingroup_set
        self.generate_group_fractions()

        # main and sub group set
        self.subgrouop_set = get_group_set(self.subgroup_list)
        self.maingroup_set = get_group_set(self.maingroup_list)

        # query r and q for group set
        self.update_group_rq()

        # query interaction between main groups
        self.query_interaction()

    def query_interaction(self):
        s_set = [str(i) for i in self.maingroup_set]
        line = ','.join(s_set)
        # connect to db file
        cwd = os.getcwd()
        db_path = os.path.join(cwd, r'data\UNIFAC.db')
        # rq_data = query_db(db_path, 'SUBGROUP, R, Q', 'GROUP_INFO', f'SUBGROUP IN ({line})')
        self.interaction_list = query_db(db_path, 'MAIN_I, MAIN_J, AIJ, AJI', 'UNIFAC_INTERACTION', f'MAIN_I IN ({line}) AND MAIN_J IN ({line})')

    def update_group_rq(self):
        # query r and q for group set
        s_set = [str(i) for i in self.subgrouop_set]
        line = ','.join(s_set)
        self.rq_dic = query_rq(line)

        for i, molecule in enumerate(self.groups_by_species):
            for j, group in enumerate(molecule['groups']):
                self.groups_by_species[i]['groups'][j]['r'] = self.rq_dic['r_dic'][group['subGroup']]
                self.groups_by_species[i]['groups'][j]['q'] = self.rq_dic['q_dic'][group['subGroup']]

    def generate_group_fractions(self):
        """
        main group fractions for each species, which is used to calculate group gamma in each species
        """
        self.group_fractions = []
        for groups in self.maingroup_list:
            fraction = np.array([group[1] for group in groups])
            self.group_fractions.append(normalize_array(fraction))

    @property
    def x(self):
        "x vector"
        return self._x

    @x.setter
    def x(self, x):
        if len(x) != len(self.groups_by_species):
            raise Exception('x list is not the same size as species list')
        return self._x


hexane = '110-54-3'
ethanol = '64-17-5'
water = '7732-18-5'
unif = UNIFAC([hexane, ethanol, water])
unif.setup_from_db()
unif.setup_r_q()
print(unif.subgroup_list)
print(unif.maingroup_list)
print(unif.maingroup_set)
print(unif.interaction_list)