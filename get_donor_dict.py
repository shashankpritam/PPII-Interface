import sys
import os
import glob
import math
import scipy
import numpy as np
import pandas as pd
import csv
import collections
from collections import defaultdict
from scipy import spatial
from scipy.spatial import KDTree

#pdb_file = sys.argv[1]

# Dictionary containing all the side chain donor atoms for all 22 amino acid residues (including ASX and GLX)
donor_dict = [('ARG', 'NE'), ('ARG', 'NH1'), ('ARG', 'NH2'), ('ASN', 'ND2'), ('ASX', 'ND2'), ('CYS', 'SG'), ('GLN', 'NE2'), ('GLX', 'NE2'), ('HIS', 'ND1'), ('HSE', 'NE2'), ('HSP', 'ND1'), ('HSP', 'NE2'), ('LYS', 'NZ'), ('SER', 'OG'), ('THR', 'OG1'), ('TRP', 'NE1'), ('TYR', 'OH')]

#print(donor_dict[0])
#print(type(donor_dict))

# Formatting the dictionary for usability
dd = collections.defaultdict(list)
for k, v in donor_dict:
    dd[k].append(v)

#for key in dd.items():
#    print(key)

#print('ASN' in dd)
#print(dd['ARG'][1])


data = []
def primary_test(pdb_file, Model_ID):
    with open(pdb_file, "r") as infile:
        for line in infile:
    #            if line.startswith("ATOM") or line.startswith("HETATM"):
            if line.startswith("ATOM"):
                # print(line.strip())
                atm_type = line[0:6].strip()
                # print(atm_type)
                atm_num = line[6:11].strip()
                atm_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                line_data = [atm_type, atm_num, atm_name, res_name, chain, res_num, x, y, z]
                data.append(line_data)
    df = pd.DataFrame(data, columns = ['atm_type', 'atm_num', 'atm_name', 'res_name', 'chain', 'res_num', 'x', 'y', 'z'])
    #print(df["x"])
    coordinates = list(zip(df['x'], df['y'], df['z']))

    atm_res_name = list(zip(df['res_name'], df['atm_name']))

    #print(coordinates)
    tree = spatial.KDTree(coordinates)
    tree.data

    data_query = df[(df['atm_name'] == 'NE1') & (df['res_name'] == 'TRP')]
    #print(data_query)

    idx = df.index[(df['atm_name'] == 'NE1') & (df['res_name'] == 'TRP')].tolist()
    #print(Model_ID, idx)

    for items in idx :
        query_result = (tree.query_ball_point(coordinates[items], r = 5.0))
        self_add = (tree.query_ball_point(coordinates[items], r = 0.0))
        neighbours = df.iloc[query_result]
        #print(neighbours)
        self_address = df.iloc[self_add]
        #print(neighbours)
        #print(neighbours.isin(self_address))
        #print(self_address)
        chain_of_trp = (self_address.at[items,'chain'])
        neighbours = neighbours.drop(neighbours[neighbours.chain != chain_of_trp].index)
        neighbours = neighbours.drop([items])
        for elements in neighbours.index:
            atm_res_name = (neighbours.res_name[elements], neighbours.atm_name[elements])
            #print(atm_res_name)
            if atm_res_name in donor_dict:
                #print(atm_res_name, neighbours.atm_num[elements], neighbours.atm_name[elements], neighbours.atm_num[elements+2], neighbours.atm_name[elements+2])
                iplus2 = (int(neighbours.atm_num[elements])+2)
                print(elements, iplus2)
                idx2 = neighbours.atm_num.index[(neighbours['atm_num'] == iplus2)].tolist()
                #print(neighbours.atm_num.index.tolist())
                #print(idx2)
                fields=[pdb_file[0:4], str(Model_ID), chain_of_trp, items, neighbours.atm_name[elements], neighbours.atm_num[elements], neighbours.atm_name.loc[elements]]
                print(fields)
                with open(r'donor_residue_info.csv', 'a') as f:
                    writer = csv.writer(f)
                    writer.writerow(fields)



#check if the file is NMR
pdb_files = glob.glob('*.pdb')
for pdb_file in pdb_files:
    output_file = pdb_file[0:4]+'_result'
    i_pdb = open(pdb_file).read().split('\n')
    i_pdb = filter(None,i_pdb)
    NMR=''
    for i in i_pdb:
        if i[:6]=='EXPDTA' and 'NMR' in i:
            NMR='true'
            break
    models=[]
    num=0
    if NMR=='true':
        i_pdb=open(pdb_file).read().split('ENDMDL')
        i_pdb=filter(None,i_pdb[:-1])
        num=0
        for i in i_pdb:
            renamed_pdb = pdb_file.split('.')[0]+'_Model_'+str(num)+'.pdb'
            w=open(renamed_pdb,'w')
            w.write(i+'\n')
            w.close()
            Model_ID = int(num)
            primary_test(renamed_pdb, Model_ID)
            num = num+1
    else:
        Model_ID = 'NA'
        primary_test(pdb_file, Model_ID)
