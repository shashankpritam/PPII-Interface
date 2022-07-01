# This program takes a pdb structure as input and return possible binding site of the PPII Helix.
# Input -
# Output -
# @Shashank Pritam - (shashankpritam@gmail.com)
# License -

# Import all the needed modules

import sys
import os
import glob
import re
from pathlib import Path
import numpy as np
import csv
import Bio.PDB
import warnings
import pandas as pd
from numpy import pi, array
from numpy.linalg import norm
from scipy.spatial.transform import Rotation as R
from Bio.PDB.vectors import rotmat
from numpy import pi
import numpy.linalg as linalg
from Bio.PDB.vectors import rotaxis2m
from Bio.PDB.vectors import Vector
from Bio.PDB.vectors import calc_angle
from Bio.PDB import PDBIO, Select
from Bio.PDB import PDBList
from Bio.PDB import PDBParser
from Bio.PDB import Selection
from Bio.PDB.vectors import Vector
from Bio.PDB.vectors import calc_angle
from Bio import BiopythonExperimentalWarning
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


current_working_dir = os.getcwd()
parser = Bio.PDB.PDBParser(QUIET=True)
pdbl = PDBList()
#input_receptor_chain_given = input("Enter the four letter Receptor PDB code along with Model ID and chain ID, eg: 1xxy0A: ")
#input_peptide_chain_given = input("Enter the four letter Peptide PDB code along with Model ID and chain ID, eg: 1xxy0B:")
# List of PDB_ID, receptor chain and peptide (ppii) chain
pdb_id = []
receptor_chain = []
ppii_chain = []
with open("pdb_chains.txt", "r") as chain_info:
    next(chain_info)
    count = 0
    for line in chain_info:
        count += 1
        pdb_id.append(line.split()[0])
        receptor_chain.append(line.split()[1])
        ppii_chain.append(line.split()[2])

#print(pdb_id, receptor_chain, ppii_chain)
#pdb_id_upper = ((input_pdb[0:4]).upper())
#pdb_id_index = pdb_id.index(pdb_id_upper)

def get_pep_chain(input_pdb):
    pdb_id_upper = (input_pdb.upper())
    pdb_id_index = pdb_id.index(pdb_id_upper)
    input_peptide_chain = ppii_chain[pdb_id_index]
    return input_peptide_chain

def get_rec_chain(input_pdb):
    pdb_id_upper = (input_pdb.upper())
    pdb_id_index = pdb_id.index(pdb_id_upper)
    input_receptor_chain = receptor_chain[pdb_id_index]
    return input_receptor_chain


# This information is needed to decide acceptor, donor, acceptor antecedent and H-Bond.
# File with list of acceptors, donors and acceptor antecedents, create dictionary for acceptor_antecedent
acceptor_file = open("H_bond_Acceptor_List").read().split('\n')
acceptor_file = filter(None,acceptor_file)
acceptor = []
acceptor_list = []
for i in acceptor_file:
    tmp=i.split('\t')
    acceptor.append(tmp[0]+':'+tmp[1])
    acceptor_list.append(tmp[-1])
# File with list of donors
donor_file = open("H_bond_Donor_List").read().split('\n')
donor_file = filter(None,donor_file)
donor = []
for i in donor_file:
    donor.append(i.replace('\t',':'))
acceptor_antecedent = {}
for atom in acceptor:
    acceptor_antecedent[atom] = acceptor_list[acceptor.index(atom)]
#print(acceptor_antecedent)
def aa_atom(atom):
    print(acceptor_list[acceptor.index(atom)])
#r_aa = acceptor_list[acceptor.index(receptor_id)]

#donor_dict is a dictionary which contains all the side donor atoms from their respective residues in a dictionary format.
donor_dict = [('ARG', 'NE'), ('ARG', 'NH1'), ('ARG', 'NH2'), ('ASN', 'ND2'), ('ASX', 'ND2'), ('CYS', 'SG'), ('GLN', 'NE2'), ('GLX', 'NE2'), ('HIS', 'ND1'), ('HSE', 'NE2'), ('HSP', 'ND1'), ('HSP', 'NE2'), ('LYS', 'NZ'), ('SER', 'OG'), ('THR', 'OG1'), ('TRP', 'NE1'), ('TYR', 'OH')]
parser = Bio.PDB.PDBParser(QUIET=True)


# Do not create a subdirectory containing PDB file
all_dataset_pdb = glob.glob('dataset/*.pdb')
hbond_files = glob.glob('data_hbond/*.txt')


# run click for each segment/pdb with all dataset pdb files
def click4all(input_pdb1, input_pdb2):
    cmd = './click '+str(input_pdb1[0])+' '+str(input_pdb2[0])+''+'>/dev/null 2>&1'
    os.system(cmd)

#Input PDB ID and get the TRP involved in Hydrogen Bond with the PPII - Only for the 39 PDB in dataset
def hbond_trp(input_pdb):
    for file in hbond_files:
        if file.split('/')[1][0:4] == input_pdb:
            list_of_trp = []
            rec_chain = get_rec_chain(input_pdb)
            with open(file, "r") as infile:
                lines = infile.readlines()
                for line in lines:
                    data = line.split(' ')
                    if "TRP" in line:
                        trp_res_id = data[data.index("TRP")+1]
                        if trp_res_id not in list_of_trp:
                            list_of_trp.append(trp_res_id)
    return list_of_trp



# This function takes input of pdb_id and an atom atom (dnr_atom, serial number = i) and its (i+2)th serial number ATOM (CX)
# and replaces that donor atom sequence id with DN and ZZ, respectively.
# suffix is used for renaming the output file.
# SCDA = Side Chain Donor ATOM
# SCDR = Side Chain Donor RESIDUE
def mask_temp_Atoms (input_pdb, scda, scdr, suffix, save_path):
    with open(input_pdb, "r") as infile:
        save_path = "click_output/"+save_path
        Path(save_path).mkdir(parents=True, exist_ok=True)
#        os.mkdir(save_path)
        file_name = input_pdb[-8:-4]+str(suffix)+'.pdb'
        completeName = os.path.join(save_path, file_name)
        the_temp_chain = get_rec_chain(input_pdb[-8:-4])
        trp_list = (hbond_trp(input_pdb[-8:-4]))
        #print(input_pdb[-8:-4], input_pdb, completeName, file_name)
        with open(completeName, 'w+') as outfile:
            lines = infile.readlines()
            for line in lines:
                if line.startswith("ATOM"):
                    atm_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain_id = line[21].strip()
                    res_seq = line[22:26].strip()
                    if res_seq in trp_list:
                        if atm_name == "NE1" and res_name == "TRP" and str(the_temp_chain) == str(chain_id):
                            line = line.replace(atm_name, "AA ", 1)
                            outfile.write(line)
                        elif atm_name == "CA" and res_name == "TRP" and str(the_temp_chain) == str(chain_id):
                            line = line.replace(atm_name, "BB", 1)
                            outfile.write(line)
                        elif atm_name == "CZ3" and res_name == "TRP" and str(the_temp_chain) == str(chain_id):
                            line = line.replace(atm_name, "CC ", 1)
                            outfile.write(line)
                        else:
                            outfile.write(line)
                    elif atm_name == scda and res_name == scdr and str(the_temp_chain) == str(chain_id):
                        if len(str(atm_name)) == 1:
                            line = line.replace(line[13:15], "NX", 1)
                            outfile.write(line)
                        elif len(str(atm_name)) == 2:
                            line = line.replace(atm_name, "NX", 1)
                            outfile.write(line)
                        elif len(str(atm_name)) == 3:
                            line = line.replace(atm_name, "NX ", 1)
                            outfile.write(line)
                        else:
                            line = line.replace(atm_name, " NX ", 1)
                            outfile.write(line)
                    elif atm_name == "CB" and res_name == scdr and str(the_temp_chain) == str(chain_id):
                        line = line.replace(atm_name, "EE", 1)
                        outfile.write(line)
                    else:
                        outfile.write(line)






def mask_query_Atoms (input_pdb, the_trp, the_nbr, the_nbr_dnr, suffix, save_path):
    save_path = "click_output/"+save_path
    Path(save_path).mkdir(parents=True, exist_ok=True)
#        os.mkdir(save_path)
    file_name = input_pdb[:-4]+str(suffix)+'.pdb'
    completeName = os.path.join(save_path, file_name)
    #print(input_pdb, completeName, file_name)
    query_structure = parser.get_structure(input_pdb[0:4], input_pdb)
    for model in query_structure:
        for chain in model:
            for residue in chain:
                if residue == the_trp:
                    the_trp_residue = residue
                    the_ne1 = the_trp_residue["NE1"]
                    the_ca = the_trp_residue["CA"]
                    the_cz3 = the_trp_residue["CZ3"]
                elif residue == the_nbr:
                    the_nbr_residue = residue
                    the_dn = the_nbr_residue[the_nbr_dnr]
                    #the_ip2 = the_nbr_residue["CX"]
    with open(input_pdb, "r") as infile:
        with open(completeName, 'w+') as outfile:
             for line in infile:
                 if line.startswith("ATOM"):
                     atm_type = line[0:6].strip()
                     atm_num = line[6:11].strip()
                     atm_name = line[12:16].strip()
                     res_name = line[17:20].strip()
                     res_seq = line[22:26].strip()
                     if atm_name == "NE1" and int(res_seq) == int(the_trp_residue.get_id()[1]):
                         line = line.replace(atm_name, "AA ", 1)
                         outfile.write(line)
                     elif atm_name == "CA" and int(res_seq) == int(the_trp_residue.get_id()[1]):
                         line = line.replace(atm_name, "BB", 1)
                         outfile.write(line)
                     elif atm_name == "CZ3" and int(res_seq) == int(the_trp_residue.get_id()[1]):
                         line = line.replace(atm_name, "CC ", 1)
                         outfile.write(line)
                     elif atm_name == the_dn.get_id() and int(res_seq) == int(the_nbr_residue.get_id()[1]):
                         if len(str(atm_name)) == 1:
                             line = line.replace(line[13:15], "NX", 1)
                             outfile.write(line)
                         elif len(str(atm_name)) == 2:
                             line = line.replace(atm_name, "NX", 1)
                             outfile.write(line)
                         elif len(str(atm_name)) == 3:
                             line = line.replace(atm_name, "NX ", 1)
                             outfile.write(line)
                         else:
                             line = line.replace(atm_name, " NX ", 1)
                             outfile.write(line)
                     elif atm_name == "CB" and int(res_seq) == int(the_nbr_residue.get_id()[1]):
                         line = line.replace(atm_name, "EE", 1)
                         outfile.write(line)
                     else:
                         outfile.write(line)
                 else:
                     outfile.write(line)


def unmask_Atoms_save (input_pdb, input_chain):
    with open(input_pdb, "r") as infile:
        output_file = input_pdb[:-4]+"_new.pdb"
        file_path = str(Path(input_pdb).resolve())
        the_dn = (file_path.split('/')[-2].split('_')[-2])
        with open(output_file, 'w+') as outfile:
             for line in infile:
                 if line.startswith("ATOM"):
                     atm_type = line[0:6].strip()
                     atm_num = line[6:11].strip()
                     atm_name = line[12:16].strip()
                     chain_id = line[21].strip()
                     if str(input_chain) == str(chain_id):
                         if atm_name == "AA":
                             line = line.replace(line[13:16], "NE1", 1)
                             outfile.write(line)
                         elif atm_name == "BB":
                             line = line.replace(line[13:15], "CA", 1)
                             outfile.write(line)
                         elif atm_name == "CC":
                             line = line.replace(line[13:16], "CZ3", 1)
                             outfile.write(line)
                         elif atm_name == "NX":
                             if len(str(the_dn)) == 1:
                                 line = line.replace(line[13:14], the_dn, 1)
                                 outfile.write(line)
                             elif len(str(the_dn)) == 2:
                                 line = line.replace(line[13:15], the_dn, 1)
                                 outfile.write(line)
                             elif len(str(the_dn)) == 3:
                                 line = line.replace(line[13:16], the_dn, 1)
                                 outfile.write(line)
                             else:
                                 line = line.replace(line[13:16], the_dn, 1)
                                 outfile.write(line)
                         elif atm_name == "EE":
                             line = line.replace(line[13:15], "CB", 1)
                             outfile.write(line)
                         else:
                             outfile.write(line)

def save_pdb (filename, pdb_id, chain):
    with open(filename, "r") as infile:
        completeName = current_working_dir+"/"+pdb_id+'_4sim.pdb'
        with open(completeName, 'w+') as outfile:
            lines = infile.readlines()
            for line in lines:
                if line.startswith("ATOM"):
                    chain_name = line[21].strip()
                    if str(chain_name) == str(chain):
                        outfile.write(line)

def change_chain (filename, chain):
        with open(filename, "r") as infile:
            pdb_id = filename[0:4]
            completeName = current_working_dir+"/"+pdb_id+"__4merge.pdb"
            with open(completeName, 'w+') as outfile:
                lines = infile.readlines()
                for line in lines:
                    if line.startswith("ATOM"):
                        chain_name = line.split()[4]
                        line = line.replace(line[21], chain)
                        outfile.write(line)


# The below function takes structure, model, chain, residue, n_atom, neighbour_atoms as input and carves that segment
# which is in the neighbourhood within 12 Angstrom of NE1 of TRP
# To carve pdb make else condition return a value of 0. It is 1 so that it selects every thing.
# Simply put - You are getting a segment of a chain from the given structure which contains TRP.
# Now each segment will have at least one another donor atom. These atoms are used for Click to find
# Binding site.
# Currently Commented Out
def carve(structure, model, chain, residue, n_atom, neighbour_atoms):
#    class AtomSelect(Select):
#        def accept_atom(self, atom):
#            if atom in neighbour_atoms:
#                return 1
#            else:
#                return 1
    pdb_id = structure.get_id()
    io = PDBIO()
    io.set_structure(structure)
# Remember to put AtomSelect() in io.save
    io.save(pdb_id+"__crvd__"+str(model.get_id())+chain.get_id()+str(residue.get_id()[1])+'_'+str(n_atom.get_parent().get_id()[1])+'.pdb')


# This functions looks for (residue, atom) pair which are within 12 Angstrom from (TRP, NE1) and
# which are not inter-hydrogen bound.
def neighbour_search(structure):
    for model in structure:
        if model.get_id() == int(0):
            for chain in model:
                for residue in chain:
                    if residue.get_resname() == 'TRP':
                        the_NE1_atom = residue["NE1"]
                        print(the_NE1_atom.get_full_id()[0], the_NE1_atom.get_full_id()[1], the_NE1_atom.get_full_id()[2], the_NE1_atom.get_full_id()[3][1], the_NE1_atom.get_parent().get_resname(), the_NE1_atom.get_full_id()[4][0], "The TRP")
                        chain_atoms  = Bio.PDB.Selection.unfold_entities(chain, 'A')
                        neighbourhood_search = Bio.PDB.NeighborSearch(chain_atoms)
                        neighbour_atoms = neighbourhood_search.search(the_NE1_atom.coord, 11.0)
    #                    print(neighbour_atoms)
    # We are once again iterating through neighbour atoms so none of the neighbour atoms are left behind.
                        for n_atom in neighbour_atoms:
                            rejection_list = []
                            atom_dic = (n_atom.get_parent().get_resname(), n_atom.get_id())
                            if (n_atom != the_NE1_atom) and (atom_dic in donor_dict):
                                n_res_atoms  = Bio.PDB.Selection.unfold_entities(n_atom.get_parent(), 'A')
# This checks if for any neighbour atom from the same chain, within 12 Angstrom except self, is side chain donor atom
# and they don't have any internal hydrogen bond (if it is then put on the rejection_list).
                                iplus2 = (n_atom.get_serial_number())+2
                                internal_look_up = neighbourhood_search.search(n_atom.coord, 3.6)
                                # We are once again already iterating through internal_atoms. This is important - from the same residue - two scda are considered!!
                                for internal_atoms in internal_look_up:
                                    if (internal_atoms != n_atom):
                                        n_atom_residue = ((n_atom.get_parent().get_resname()))
                                        internal_look_up_residue = ((internal_atoms.get_parent().get_resname()))
                                        n_atom_id = (n_atom.get_name()+":"+n_atom_residue)
                                        internal_look_up_residue_id = (internal_atoms.get_name()+":"+internal_look_up_residue)
                                        if (n_atom_id in donor and internal_look_up_residue_id in acceptor):
                                            i_aa_atom = internal_atoms.get_parent()[acceptor_antecedent[internal_look_up_residue_id]]
                                            n_vector = n_atom.get_vector()
                                            i_vector = internal_atoms.get_vector()
                                            i_aa_vector = i_aa_atom.get_vector()
                                            i_angle_hbond = (calc_angle(n_vector, i_vector, i_aa_vector))
                                            if np.degrees(i_angle_hbond) >90.0 and np.degrees(i_angle_hbond)<180:
                                                rejection_list.append(n_atom)
                                                #print(np.degrees(i_angle_hbond))
                                        elif (internal_look_up_residue_id in donor and n_atom_id in acceptor):
                                            n_aa_atom = n_atom.get_parent()[acceptor_antecedent[n_atom_id]]
                                            n_vector = n_atom.get_vector()
                                            i_vector = internal_atoms.get_vector()
                                            n_aa_vector = n_aa_atom.get_vector()
                                            n_angle_hbond = (calc_angle(i_vector, n_vector, n_aa_vector))
                                            if np.degrees(n_angle_hbond) >90.0 and np.degrees(n_angle_hbond)<180:
                                                rejection_list.append(n_atom)
                                                #print(np.degrees(n_angle_hbond))
                                        else:
                                            continue
                                        #print(rejection_list) #The internal hydrogen bond look up ends here
                                if (n_atom not in rejection_list) and (n_atom != the_NE1_atom) and (atom_dic in donor_dict):
                                    #Change this if you need to. One representative atom from another side chain donor
                                    CX = str("CB")
    #                                        for na in n_res_atoms:
    #                                            if (na.get_serial_number()) == iplus2:
    #                                                CX = na.get_id()
    #                                            else:
    #                                                print("something is not right")
                                    representative_atom = (n_atom.get_id())
                                    #print(n_atom.get_full_id()[0], n_atom.get_full_id()[1], n_atom.get_full_id()[2], n_atom.get_full_id()[3][1], n_atom.get_parent().get_resname(), n_atom.get_full_id()[4][0], the_NE1_atom-n_atom)
                                    carve(structure, model, chain, residue, n_atom, neighbour_atoms)
                                    for dataset_file in all_dataset_pdb:
                                        # Very important - format of pdb file is dataset/XXYY.pdb, hence int(16)
                                        if len(dataset_file) == int(16): #and (dataset_file == "dataset/1gbq.pdb"):
                                            #print(dataset_file)
                                            the_path = str(n_atom.get_full_id()[0])+"_"+str(n_atom.get_full_id()[1])+"_"+str(n_atom.get_full_id()[2])+"_"+str(the_NE1_atom.get_full_id()[3][1])+"_"+str(n_atom.get_full_id()[3][1])+"_"+str(n_atom.get_parent().get_resname()+"_"+n_atom.get_id()+"_"+str(dataset_file[-8:-4]))
                                            #print(the_path)
                                            mask_temp_Atoms(dataset_file, str(n_atom.get_id()), str(n_atom.get_parent().get_resname()), "_rnmd_ds", the_path)
                                            mask_query_Atoms (structure.get_id()+"__crvd__"+str(model.get_id())+chain.get_id()+str(residue.get_id()[1])+'_'+str(n_atom.get_parent().get_id()[1])+'.pdb', residue, n_atom.get_parent(), n_atom.get_id(), '__rnmd', the_path)
                                            #mask_Atoms(structure.get_id()+"__crvd__"+str(model.get_id())+chain.get_id()+str(residue.get_id()[1])+'_'+str(n_atom.get_parent().get_id()[1])+'.pdb', representative_atom,  CX, '__rnmd', the_path)

                                #return representative_atom





#Input Output and Command
#init_ds = open("pdb_44_dataset.txt").read().split('\n')
#for pdb in init_ds:
#input_pdb_given = pdb
input_pdb_given = input("Enter the four letter PDB code of Query Protein: ")
pdbl = PDBList()
print("Input PDB Given : "+input_pdb_given)
input_pdb = pdbl.retrieve_pdb_file(input_pdb_given, obsolete=False, pdir=None, file_format="pdb", overwrite=False)
input_structure = parser.get_structure(input_pdb_given, input_pdb)
neighbour_search(input_structure)


files_4_click = glob.glob(current_working_dir+"/click_output/*/", recursive = True)
for folders in files_4_click:
    if folders.startswith(current_working_dir+"/click_output/"+input_pdb_given.upper()) or folders.startswith(current_working_dir+"/click_output/"+input_pdb_given.lower()):
        #print(folders)
        renamed_pdb = glob.glob(folders+'/*_rnmd.pdb')
        dataset_renamed_file = glob.glob(folders+'/*_rnmd_ds.pdb')
        click4all(renamed_pdb, dataset_renamed_file)

click_atoms = ["AA", "BB", "CC", "NX", "EE"]
predicted_alignments = []
for folders in files_4_click:
    if folders.startswith(current_working_dir+"/click_output/"+input_pdb_given.upper()) or folders.startswith(current_working_dir+"/click_output/"+input_pdb_given.lower()):
        click_file = glob.glob(folders+'/*.clique')
        #print(click_file)
        with open(click_file[0], 'r') as infile:
            data_to_read = infile.read()
            data = data_to_read.split()
            split_data = data_to_read.split()
            Matched_Atoms = int(split_data[6])
            RMSD = float(split_data[9])
            Structure_Overlap = float(split_data[13])
            #print(RMSD, Matched_Atoms)
            #print(click_file[0])
            carved_frag_info = (click_file[0]).split('/')[-1]
            #print(carved_frag_info, RMSD, Structure_Overlap, Matched_Atoms)
            if RMSD <= 0.6 and Matched_Atoms == 4:
                predicted_alignments.append(carved_frag_info+"_"+str(RMSD)+"_"+str(Structure_Overlap))
                #print(carved_frag_info, RMSD, Structure_Overlap, Matched_Atoms, "First Criteria")
            elif RMSD <= 1.0 and Matched_Atoms == 5:
                predicted_alignments.append(carved_frag_info+"_"+str(RMSD)+"_"+str(Structure_Overlap))
                #print(carved_frag_info, RMSD, Structure_Overlap, Matched_Atoms, "Second Criteria")

#print(predicted_alignments[0].split('_'))
list_of_unique_alignment = []
for alignments in predicted_alignments:
    #print(alignments)
    alignments = alignments.split('_')
    #print(alignments, alignments[-1], alignments[-2])
    #print(alignments[0], (alignments[7][-4:]))
    if (alignments[0]).casefold() != (alignments[7][-4:]).casefold():
        list_of_unique_alignment.append(alignments)

best_alignment = sorted(list_of_unique_alignment, key = lambda x: (-float(x[-2]), -float(x[-1])))[0]
#for item in list_of_unique_alignment:
#    print(item)

print("For query structure ", best_alignment[0], ", predicted binding site is at - ", "Model = ", best_alignment[4][0],  "Chain = ",best_alignment[4][1], "TRP = ", best_alignment[4][2:], "NBR = ", best_alignment[5])
print("Best alignment details are - ", "PDB ID = ", best_alignment[7][-4:],  "RMSD = ", best_alignment[-2], "SO = ", best_alignment[-1])




def run_sim_sim (input_pdb):
    for folders in files_4_click:
        if folders.startswith(current_working_dir+"/click_output/"+input_pdb_given.upper()) or folders.startswith(current_working_dir+"/click_output/"+input_pdb_given.lower()):
            dataset_renamed_file = glob.glob(folders+'/*_rnmd.1.pdb')
            renamed_pdb = glob.glob(folders+'/*_rnmd_ds.1.pdb')
            click_file = glob.glob(folders+'/*.clique')
            carved_frag_info = (click_file[0]).split('/')[8]
            if carved_frag_info.split("_") == best_alignment[:-2]:
                #print(best_alignment, click_file[0])
                receptor_chain = best_alignment[4][1]
                peptide_chain = get_pep_chain(best_alignment[7][-4:])
                #print("receptor chain = ", receptor_chain, "peptide chain = ", peptide_chain)
                unmask_Atoms_save(renamed_pdb[0], receptor_chain)
                unmask_Atoms_save(dataset_renamed_file[0], peptide_chain)
                renamed_pdb_unmask = renamed_pdb[0][:-4]+"_new.pdb"
                dataset_renamed_file_unmask = dataset_renamed_file[0][:-4]+"_new.pdb"

                save_pdb (renamed_pdb_unmask, best_alignment[0], receptor_chain)
                save_pdb (dataset_renamed_file_unmask, best_alignment[7][-4:], peptide_chain)

                predicted_receptor_structure = parser.get_structure(best_alignment[0], renamed_pdb_unmask)
                template_peptide_structure = parser.get_structure((best_alignment[7][-4:]), dataset_renamed_file_unmask)

                input_receptor_chain_given = predicted_receptor_structure.get_id()+str(0)+receptor_chain
                input_peptide_chain_given = template_peptide_structure.get_id()+str(0)+peptide_chain
    return (predicted_receptor_structure, receptor_chain, template_peptide_structure, peptide_chain)

#print(run_sim_sim(input_pdb_given))


parser = Bio.PDB.PDBParser(QUIET=True)
pdbl = PDBList()

input_receptor_chain_given = run_sim_sim(input_pdb_given)[0].get_id()
input_peptide_chain_given = run_sim_sim(input_pdb_given)[2].get_id()

input_receptor_chain = run_sim_sim(input_pdb_given)[1]
input_peptide_chain = run_sim_sim(input_pdb_given)[3]

#Provide Model of relevance here
input_receptor_model = 0



input_receptor_structure = run_sim_sim(input_pdb_given)[0]
input_peptide_structure = run_sim_sim(input_pdb_given)[2]



print("Receptor PDB ID : "+input_receptor_chain_given)
print("Peptide PDB ID : "+input_peptide_chain_given)
print(input_receptor_chain_given, input_receptor_chain, input_peptide_chain_given, input_peptide_chain)


class ReceptorSelect(Select):
    def accept_model(self, model):
        if model.get_id() == int(input_receptor_model):
            return 1
        else:
            return 0
    def accept_residue(self, residue):
        if residue.id[0] == " ":
            return 1
        else:
            return 0

class PeptideSelect(Select):
    def accept_residue(self, residue):
        if residue.id[0] == " ":
            return 1
        else:
            return 0

Main_Chain_Atoms = ["C", "N", "O", "CA"]
class MCAtomSelect(Select):
    def accept_atom(self, atom):
        if atom.get_id() in Main_Chain_Atoms:
            return 1
        else:
            return 0

io = PDBIO()
#io.set_structure(input_receptor_structure)
#io.save(input_receptor_chain_given+"_receptor.pdb")
#io.set_structure(input_peptide_structure)
#io.save(input_peptide_chain_given+"_peptide.pdb")

rec_structure = parser.get_structure(input_receptor_chain_given, input_receptor_chain_given+"_4sim.pdb")
pep_structure = parser.get_structure(input_peptide_chain_given, input_peptide_chain_given+"_4sim.pdb")

print(rec_structure.get_id())
'''for model_rec in rec_structure:
    for chain_rec in model_rec:
        chain_rec.id = 'X'
        io.set_structure(rec_structure)
        io.save(input_receptor_chain_given+"__4merge.pdb")

for model_pep in pep_structure:
    for chain_pep in model_pep:
        chain_pep.id = 'Y'
        io.set_structure(pep_structure)
        io.save(input_peptide_chain_given+"__4merge.pdb")'''


change_chain(input_receptor_chain_given+"_4sim.pdb" , 'X')
change_chain(input_peptide_chain_given+"_4sim.pdb" , 'Y')
new_rec_structure = parser.get_structure(input_receptor_chain_given, input_receptor_chain_given+"__4merge.pdb")
new_pep_structure = parser.get_structure(input_peptide_chain_given, input_peptide_chain_given+"__4merge.pdb")

rec_file = input_receptor_chain_given+"__4merge.pdb"
pep_file = input_peptide_chain_given+"__4merge.pdb"

combined_file = input_receptor_chain_given+"__"+input_peptide_chain_given+"__4sim.pdb"
#print(rec_file, pep_file, combined_file)


with open(combined_file, 'w') as outfile:
    with open(rec_file) as rec_file:
        for line in rec_file:
            if not line.startswith('END'):
                outfile.write(line)
    with open(pep_file) as pep_file:
            for line in pep_file:
                outfile.write(line)

command4 = input_receptor_chain_given+"__4merge.pdb"
command5 = input_peptide_chain_given+"__4merge.pdb"
command7 = input_receptor_chain_given+"_receptor.pdb"
command8 = input_peptide_chain_given+"_peptide.pdb"
#os.remove(command4)
#os.remove(command5)
#os.remove(command7)
#os.remove(command8)

simulation_pdb = parser.get_structure(input_receptor_chain_given[0:4], input_receptor_chain_given+"__"+input_peptide_chain_given+"__4sim.pdb")


#Hydrogen bond distance
distance = 3.5
Main_Chain_Atoms = ["C", "N", "O", "CA"]
def NHBS(structure):
    for model in structure:
        id = model.get_id()
        for chain in model:
            if chain.get_id() == "X":
                r_chain = chain
            elif chain.get_id() == "Y":
                p_chain = chain
            else:
                print("something is not right")
        receptor_atoms  = Bio.PDB.Selection.unfold_entities(r_chain, 'A')
        peptide_atoms = Bio.PDB.Selection.unfold_entities(p_chain, 'A')
        hb = int(0)
        #print(donor, acceptor)
        for p_atoms in peptide_atoms:
            for r_atoms in receptor_atoms:
                if (p_atoms)-(r_atoms) <= distance:
                    #print ((p_atoms.get_full_id()))
                    peptide_residue = ((p_atoms.get_parent().get_resname()))
                    receptor_residue = ((r_atoms.get_parent().get_resname()))
                    peptide_id = (p_atoms.get_name()+":"+peptide_residue)
                    receptor_id = (r_atoms.get_name()+":"+receptor_residue)
                    #print(peptide_id, receptor_id)
                    #print(p_atoms.get_vector())
                    if (peptide_id in donor and receptor_id in acceptor):
                        r_aa_atom = r_atoms.get_parent()[acceptor_antecedent[receptor_id]]
                        p_vector = p_atoms.get_vector()
                        r_vector = r_atoms.get_vector()
                        r_aa_vector = r_aa_atom.get_vector()
                        r_angle_hbond = (calc_angle(p_vector, r_vector, r_aa_vector))
                        if np.degrees(r_angle_hbond) >100.0 and np.degrees(r_angle_hbond)<180:
                            hb = hb+1
                    elif (receptor_id in donor and peptide_id in acceptor):
                        p_aa_atom = p_atoms.get_parent()[acceptor_antecedent[peptide_id]]
                        p_vector = p_atoms.get_vector()
                        r_vector = r_atoms.get_vector()
                        p_aa_vector = p_aa_atom.get_vector()
                        p_angle_hbond = (calc_angle(p_aa_vector, p_vector, r_vector))
                        if np.degrees(p_angle_hbond) >100.0 and np.degrees(p_angle_hbond)<180:
                            hb = hb+1
                    else:
                        hb = hb+0
        #print((hb**2))
        return (hb**2)

#print(NHBS(simulation_pdb))

lower_bound_dist = float(2.8)
upper_bound_dist = float(4.2)
lower_bound_ang = (80)
upper_bound_ang = (180)
clash_dist = float(2.8)
PVWS_upper = float(5.0)
temperature = float(400.0)#float(673.15) In Kelvin
gas_constant = float(8.31) #In J⋅K−1⋅mol−1
translation = float(2) #± 2 Å
rotation = float(5) #± 5°
scale_distance_factor = float(1)

def RHBS(structure):
    for model in structure:
        id = model.get_id()
        for chain in model:
            if chain.get_id() == "X":
                r_chain = chain
            elif chain.get_id() == "Y":
                p_chain = chain
            else:
                print("something is not right")
        receptor_atoms  = Bio.PDB.Selection.unfold_entities(r_chain, 'A')
        peptide_atoms = Bio.PDB.Selection.unfold_entities(p_chain, 'A')
        dist_score = 0
        ang_score = 0
        #print(donor, acceptor)
        for p_atoms in peptide_atoms:
            for r_atoms in receptor_atoms:
                peptide_residue = ((p_atoms.get_parent().get_resname()))
                receptor_residue = ((r_atoms.get_parent().get_resname()))
                peptide_id = (p_atoms.get_name()+":"+peptide_residue)
                receptor_id = (r_atoms.get_name()+":"+receptor_residue)
                if (peptide_id in donor and receptor_id in acceptor):
                    dist = p_atoms-r_atoms
                    r_aa_atom = r_atoms.get_parent()[acceptor_antecedent[receptor_id]]
                    p_vector = p_atoms.get_vector()
                    r_vector = r_atoms.get_vector()
                    r_aa_vector = r_aa_atom.get_vector()
                    r_angle_hbond = (calc_angle(p_vector, r_vector, r_aa_vector))
                    if (dist > lower_bound_dist) and (dist < upper_bound_dist):
                        #print(dist)
                        dist_score = dist_score + scale_distance_factor*(dist-lower_bound_dist)*(dist-upper_bound_dist)
                        #print (dist_score)
                    else:
                        dist_score = dist_score + 0
                    if np.degrees(r_angle_hbond)>lower_bound_ang and np.degrees(r_angle_hbond)<upper_bound_ang:
                        #print(np.degrees(r_angle_hbond))
                        ang_score = ang_score + (np.degrees(r_angle_hbond)-lower_bound_ang)*(np.degrees(r_angle_hbond)-upper_bound_ang)
                        ang_score = np.radians(ang_score)
                        #print(ang_score)
                    else:
                        ang_score = ang_score + 0
                elif (receptor_id in donor and peptide_id in acceptor):
                    dist = p_atoms-r_atoms
                    p_aa_atom = p_atoms.get_parent()[acceptor_antecedent[peptide_id]]
                    p_vector = p_atoms.get_vector()
                    r_vector = r_atoms.get_vector()
                    p_aa_vector = p_aa_atom.get_vector()
                    p_angle_hbond = (calc_angle(p_aa_vector, p_vector, r_vector))
                    if (dist > lower_bound_dist) and (dist < upper_bound_dist):
                        dist_score = dist_score + scale_distance_factor*(dist-lower_bound_dist)*(dist-upper_bound_dist)
                        #print (dist_score)
                    else:
                        dist_score = dist_score + 0
                    if np.degrees(p_angle_hbond)>lower_bound_ang and np.degrees(p_angle_hbond)<upper_bound_ang:
                        ang_score = ang_score + (np.degrees(p_angle_hbond)-lower_bound_ang)*(np.degrees(p_angle_hbond)-upper_bound_ang)
                        ang_score = np.radians(ang_score)
                        #print(ang_score)
                    else:
                        ang_score = ang_score + 0
                else:
                    continue
        return (dist_score+(ang_score))
        #return (dist_score+(4.6*(ang_score)))

#RHBS(simulation_pdb)

def CS(structure):
    for model in structure:
        id = model.get_id()
        for chain in model:
            if chain.get_id() == "X":
                r_chain = chain
            elif chain.get_id() == "Y":
                p_chain = chain
            else:
                print("something is not right")
        receptor_atoms  = Bio.PDB.Selection.unfold_entities(r_chain, 'A')
        peptide_atoms = Bio.PDB.Selection.unfold_entities(p_chain, 'A')
#        clash_score = 0
        MCMC = 0
        MCSC = 0
        SCSC = 0
        #print(donor, acceptor)
        for p_atoms in peptide_atoms:
            for r_atoms in receptor_atoms:
                dist = p_atoms-r_atoms
#                if dist < clash_dist:
#                    print(dist, p_atoms, r_atoms)
#                    dist_value = dist**2
#                    dist_value = (dist_value+1)**2
#                else:
#                    dist_value = dist_value+0
                if (dist < clash_dist) and (r_atoms.get_id() in Main_Chain_Atoms) and (p_atoms.get_id() in Main_Chain_Atoms):
                    MCMC = MCMC+1
#                    clash_score = clash_score + (100*dist_value)
                elif (dist < clash_dist) and (r_atoms.get_id() not in Main_Chain_Atoms) and (p_atoms.get_id() in Main_Chain_Atoms):
                    MCSC = MCSC+1
#                    clash_score = clash_score + (5*((dist_value**2)**0.5))
                elif (dist < clash_dist) and (r_atoms.get_id() in Main_Chain_Atoms) and (p_atoms.get_id() not in Main_Chain_Atoms):
                    MCSC = MCSC+1
#                    clash_score = clash_score + (5*((dist_value**2)**0.5))
                elif (dist < clash_dist) and (r_atoms.get_id() not in Main_Chain_Atoms) and (p_atoms.get_id() not in Main_Chain_Atoms):
                    SCSC = SCSC+1
#                    clash_score = clash_score + (5*(dist_value))
                else:
                    continue
        #clash_score = (100*(MCMC**2))+(5*(MCSC**2))+(5*(SCSC))
        clash_score = (10*(MCMC))+((MCSC**0.5))+((SCSC*0.25))
        #print(MCMC, MCSC, SCSC, clash_score)
        return (clash_score)

#CS(simulation_pdb)

def PVWS(structure):
    for model in structure:
        id = model.get_id()
        for chain in model:
            if chain.get_id() == "X":
                r_chain = chain
            elif chain.get_id() == "Y":
                p_chain = chain
            else:
                print("something is not right")
        receptor_atoms  = Bio.PDB.Selection.unfold_entities(r_chain, 'A')
        peptide_atoms = Bio.PDB.Selection.unfold_entities(p_chain, 'A')
        PVWS_num = 0
        #print(donor, acceptor)
        for p_atoms in peptide_atoms:
            for r_atoms in receptor_atoms:
                dist = p_atoms-r_atoms
                if (dist > clash_dist) and (dist < PVWS_upper):
                    #print(dist)
                    PVWS_num = PVWS_num+1
                else:
                    PVWS_num = PVWS_num+0
        return(PVWS_num)
        #return (PVWS_num)

#PVWS(simulation_pdb)
def total_energy(simulation_pdb):
    energy = -(NHBS(simulation_pdb))+(RHBS(simulation_pdb))+(CS(simulation_pdb))-(PVWS(simulation_pdb))
    return(energy)
#print("Total Energy =", total_energy(simulation_pdb), ", NHBS = ", (NHBS(simulation_pdb)), ", RHBS = ", (RHBS(simulation_pdb)), ", CS = ", (CS(simulation_pdb)), ", PVWS = ", (PVWS(simulation_pdb)))

'''NHBS(simulation_pdb)
RHBS(simulation_pdb)
CS(simulation_pdb)
PVWS(simulation_pdb)
'''

def metropolis_p(E1,E2):
    diff = float(E1-E2)
    div = float(gas_constant*temperature)
    diff_by_div = diff/div
    P = np.exp(diff_by_div)
    return P

#print(metropolis_p(2700000, total_energy))
def x_trans(atom, dx, rot):
    return atom.coord+(np.array([float(dx), 0, 0], dtype=float))
    #return np.dot(atom.coord, rot)+(np.array([float(dx), 0, 0], dtype=float))
def y_trans(atom, dy, rot):
    return atom.coord+(np.array([0, float(dy), 0], dtype=float))
    #return np.dot(atom.coord, rot)+(np.array([0, float(dy), 0], dtype=float))
def z_trans(atom, dz, rot):
    return atom.coord+(np.array([0, 0, float(dz)], dtype=float))
    #return np.dot(atom.coord, rot)+(np.array([0, 0, float(dz)], dtype=float))

def calc_angle_vectors(vector1, vector2):
    vec_1 = vector1/np.linalg.norm(vector1)
    vec_2 = vector2/np.linalg.norm(vector2)
    dot_product = np.dot(vec_1, vec_2)
    angle = np.arccos(dot_product)
    return np.degrees(angle)

def rot_one_on_two(vector1, vector2, angle):
    comp_parallel = vector2*(np.dot(vector1, vector2) / np.dot(vector2, vector2))
    comp_normal = vector1 - comp_parallel
    b_cross_a_norm = np.cross(vector2, comp_normal)
    rotated_vector = comp_parallel + (comp_normal * np.cos(angle))+(linalg.norm(comp_normal) * (b_cross_a_norm / linalg.norm(b_cross_a_norm)) * np.sin(angle))
    return (rotated_vector)

def move_pdb(structure, iteration):
    for model in structure:
        for chain in model:
            if chain.get_id() == "Y":
                the_end_list = np.empty((0,3), float)
                the_peptide_atoms = Bio.PDB.Selection.unfold_entities(chain, 'A')
                for atom in the_peptide_atoms:
                    atom_coord = np.array(atom.coord)
                    the_end_list = np.append(the_end_list, [atom_coord], axis=0)
                #print(the_end_list)
                the_first_end = the_end_list[0]
                the_last_end = the_end_list[len(the_end_list)-1]
                # mean_coordinate is the mean coordinate through all the peptide atoms
                mean_coordinate = (np.sum(the_end_list, axis=0)/len(the_end_list))
                direction_vector = (the_first_end-the_last_end)
                direction_vector = direction_vector/ np.sqrt(np.sum(direction_vector**2))
                random_vector = mean_coordinate #+ ((np.random.uniform(-1, 1))*direction_vector)
                random_vector = random_vector/norm(random_vector)
                rotation_degrees = np.radians(45)
                a = np.array([0,1,0])
                b = np.array([1,0,0])
                rotated_vector = rot_one_on_two(random_vector, a, rotation_degrees)
                #print(a, b, random_vector, rotated_vector)
                #print(calc_angle_vectors(a,b), calc_angle_vectors(a, rotated_vector), calc_angle_vectors(b, rotated_vector), calc_angle_vectors(random_vector, rotated_vector))
                dx = 2*(np.random.uniform(-0.5, 0.5))
                dy = 2*(np.random.uniform(-0.5, 0.5))
                dz = 2*(np.random.uniform(-0.5, 0.5))
                rot = 2*(np.random.uniform(-0.5, 0.5, (3, 3)))
                #print(rot)
                #print(dx, dy, dz)
                for residue in chain:
                    for atom in residue:
                        #print(dx, dy, dz, rot)
                        atom.coord = (x_trans(atom, dx, rot))
                        atom.coord = (y_trans(atom, dy, rot))
                        atom.coord = (z_trans(atom, dz, rot))
                        #print(x_trans(atom, dx, rot), y_trans(atom, dy, rot), z_trans(atom, dz, rot))
        #io = PDBIO()
        io.set_structure(structure)
        file_name = "pdb_traj"+str(iteration)+".pdb"
        io.save(file_name)

energy_data = [0]
energy_plot = []
iteration_axis = []
#To move the coordinates of the structure
def decide_move(input_structure, iteration):
    limit = int(20)
    E2_Val = float(total_energy(input_structure))
    last_energy_position = (len(energy_data)-1)
    E1_Val = float(energy_data[last_energy_position])
    while iteration <= limit:
        if iteration == 0:
            energy_data[0] = (E2_Val)
            #print(E2_Val)
            iteration = iteration+1
            move_pdb(input_structure, iteration)
        else:
            traj_file_name = "pdb_traj"+str(iteration)+".pdb"
            traj_file_str = parser.get_structure(input_receptor_chain_given[0:4], traj_file_name)
            saved_file_name = "pdb_traj_saved"+str(iteration)+".pdb"
            P_Val = (metropolis_p(E1_Val, E2_Val))
            #print(P_Val)
            if P_Val > 1:
                print("Len8OPS", P_Val, E2_Val, iteration)
                energy_data.append(E2_Val)
                io.set_structure(input_structure)
                io.save(saved_file_name)
#                print("Len8OPS", P_Val, E2_Val, iteration)
                energy_plot.append(E2_Val)
                iteration_axis.append(iteration)
                iteration = iteration+1
                move_pdb(traj_file_str, iteration)
                decide_move(traj_file_str, iteration)
            else:
                print("Not_OPS", P_Val, E2_Val, iteration)
                iteration = iteration+0
                move_pdb(input_structure, iteration)
                decide_move(input_structure, iteration)
            break


#decide_move(simulation_pdb, 0)

##To create a single multi-model .pdb file - out_combined.pdb
with open("out_combined.pdb", 'w+') as outfile:
    for files in glob.glob('pdb_traj_saved?.pdb'):
        with open(files) as file:
                for line in file:
                    outfile.write(line)
command12 = "sed -i 's/END/ENDMDL/g' out_combined.pdb"
#os.system(command12)


#Important Snippet
