import sys
import os
import glob
import pandas as pd
import csv
from scipy.spatial import KDTree

# Dictionary containing all the side chain donor atoms for all 22 amino acid residues (including ASX and GLX)
donor_dict = [('ARG', 'NE'), ('ARG', 'NH1'), ('ARG', 'NH2'), ('ASN', 'ND2'), ('ASX', 'ND2'), ('CYS', 'SG'),
              ('GLN', 'NE2'), ('GLX', 'NE2'), ('HIS', 'ND1'), ('HSE', 'NE2'), ('HSP', 'ND1'), ('HSP', 'NE2'),
              ('LYS', 'NZ'), ('SER', 'OG'), ('THR', 'OG1'), ('TRP', 'NE1'), ('TYR', 'OH')]

def primary_test(pdb_file, Model_ID):
    data = []
    with open(pdb_file, "r") as infile:
        for line in infile:
            if line.startswith("ATOM"):
                atm_type = line[0:6].strip()
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

    df = pd.DataFrame(data, columns=['atm_type', 'atm_num', 'atm_name', 'res_name', 'chain', 'res_num', 'x', 'y', 'z'])
    coordinates = df[['x', 'y', 'z']].values

    tree = KDTree(coordinates)

    data_query = df[(df['atm_name'] == 'NE1') & (df['res_name'] == 'TRP')]
    idx = data_query.index.tolist()

    with open('donor_residue_info.csv', 'a') as f:
        writer = csv.writer(f)
        for idx in idx:
            query_result = tree.query_ball_point(coordinates[idx], r=5.0)
            self_add = tree.query_ball_point(coordinates[idx], r=0.0)
            neighbours = df.iloc[query_result]
            self_address = df.iloc[self_add]
            chain_of_trp = self_address.at[idx, 'chain']
            neighbours = neighbours[neighbours['chain'] == chain_of_trp].drop(idx)
            for elements in neighbours.index:
                atm_res_name = (neighbours.at[elements, 'res_name'], neighbours.at[elements, 'atm_name'])
                if atm_res_name in donor_dict:
                    iplus2 = int(neighbours.at[elements, 'atm_num']) + 2
                    idx2 = neighbours[neighbours['atm_num'] == str(iplus2)].index.tolist()
                    fields = [pdb_file[0:4], str(Model_ID), chain_of_trp, idx, neighbours.atm_name[elements], neighbours.atm_num[elements], neighbours.atm_name.loc[elements]]
                    writer.writerow(fields)

# Check if the file is NMR
pdb_files = glob.glob('*.pdb')
for pdb_file in pdb_files:
    output_file = pdb_file[0:4] + '_result'
    i_pdb = open(pdb_file).read().split('\n')
    i_pdb = filter(None, i_pdb)
    NMR = any('NMR' in i for i in i_pdb if i[:6] == 'EXPDTA')
    if NMR:
        i_pdb = open(pdb_file).read().split('ENDMDL')
        i_pdb = filter(None, i_pdb[:-1])
        for num, i in enumerate(i_pdb):
            renamed_pdb = pdb_file.split('.')[0] + '_Model_' + str(num) + '.pdb'
            with open(renamed_pdb, 'w') as w:
                w.write(i + '\n')
            Model_ID = num
            primary_test(renamed_pdb, Model_ID)
    else:
        Model_ID = 'NA'
        primary_test(pdb_file, Model_ID)
