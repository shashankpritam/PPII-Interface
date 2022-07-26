import Bio.PDB
import glob
import numpy as np
from pathlib import Path
from Bio.PDB.vectors import calc_angle
from Bio.PDB import PDBIO, Select, PDBList
from modules.carve_pdb import carve
from modules.mask_temp_Atoms import mask_temp_Atoms
from modules.mask_query_Atoms import mask_query_Atoms
# Subdirectories containing the template PDB files and the hydrogen bond info.
all_dataset_pdb = glob.glob('dataset/*.pdb')
# This hydrogen bond info had 4 Angstrom as distance cut-off between the acceptor and the donor atoms.
# The range for hydrogen bond angle was (90, 180).
hbond_files = glob.glob('data_hbond/*.txt')
# This information is needed to decide acceptor, donor, acceptor antecedent atoms for hydrogen bond inference.
# Loads file with list of acceptors, donors and acceptor antecedents, create dictionary for acceptor_antecedent

# File with list of acceptors
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

# Acceptor Antecedent Dictionary
acceptor_antecedent = {}
for atom in acceptor:
    acceptor_antecedent[atom] = acceptor_list[acceptor.index(atom)]

# A function to get the acceptor antecedent atom of an acceptor atom involved in hydrogen bond.
def aa_atom(atom):
    print(acceptor_list[acceptor.index(atom)])
# This functions below - neighbour_search looks for a side chain donor atom for hydrogen bond which are within 12 Angstrom from any (TRP, NE1)
# and considers only those which are not inter-hydrogen bound within the same chain.
# Please do not set/change the parameters here.
# Input is simply a biopython structure Object
# Neighbourhood Search Parameters
parameter_file = open("param.txt").read().split('\n')
neighbourhood_look_up_cut_off = float(parameter_file[5].split( )[2])
H_Bond_Cut_Off = float(parameter_file[6].split( )[2])
# donor_dict is a dictionary which contains all the side donor atoms from their respective residues in a dictionary format.
# The information of donor_dict is taken from the Modeller software.
donor_dict = [('ARG', 'NE'), ('ARG', 'NH1'), ('ARG', 'NH2'), ('ASN', 'ND2'), ('ASX', 'ND2'), ('CYS', 'SG'), ('GLN', 'NE2'), ('GLX', 'NE2'), ('HIS', 'ND1'), ('HSE', 'NE2'), ('HSP', 'ND1'), ('HSP', 'NE2'), ('LYS', 'NZ'), ('SER', 'OG'), ('THR', 'OG1'), ('TRP', 'NE1'), ('TYR', 'OH')]


def neighbour_search(structure):
    for model in structure:
        if model.get_id() == int(0):
            for chain in model:
                for residue in chain:
# Looks if in the given structure if TRP is present on any postion (iterating through model and chain)
# If you want to provide Model ID for your query (NMR structure) change the line model.get_id() == int(0):
# as per your requirements. Else deault model is set to 0.
                    if residue.get_resname() == 'TRP':
                        the_NE1_atom = residue["NE1"]
# The TRP info will be shown as console output.
                        print("Residue Tryptophan is present at : ", the_NE1_atom.get_full_id()[0], the_NE1_atom.get_full_id()[1], the_NE1_atom.get_full_id()[2], the_NE1_atom.get_full_id()[3][1])
                        model_atoms  = Bio.PDB.Selection.unfold_entities(model, 'A')
                        neighbourhood_search = Bio.PDB.NeighborSearch(model_atoms)
                        neighbour_atoms = neighbourhood_search.search(the_NE1_atom.coord, neighbourhood_look_up_cut_off)
# Saves all the neighbour atoms within cut off of neighbourhood_look_up_cut_off
# Now we are iterating through all the neighbour atoms so none of the neighbour atoms are left behind.
                        for n_atom in neighbour_atoms:
                            rejection_list = []
                            atom_dic = (n_atom.get_parent().get_resname(), n_atom.get_id())
                            if (n_atom != the_NE1_atom) and (atom_dic in donor_dict):
                                #print(the_NE1_atom.get_full_id(), n_atom.get_full_id())
#                                n_res_atoms  = Bio.PDB.Selection.unfold_entities(n_atom.get_parent(), 'A')
# This checks if for any neighbour atom from the same chain, within 12 Angstrom except self, is side chain donor atom
# and they don't have any internal hydrogen bond (if it is then put on the rejection_list).
#                                iplus2 = (n_atom.get_serial_number())+2
                                internal_look_up = neighbourhood_search.search(n_atom.coord, H_Bond_Cut_Off)
# We are once again already iterating through internal_atoms. This is important - from the same residue - two scda might be considered!!
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
                                        elif (internal_look_up_residue_id in donor and n_atom_id in acceptor):
                                            n_aa_atom = n_atom.get_parent()[acceptor_antecedent[n_atom_id]]
                                            n_vector = n_atom.get_vector()
                                            i_vector = internal_atoms.get_vector()
                                            n_aa_vector = n_aa_atom.get_vector()
                                            n_angle_hbond = (calc_angle(i_vector, n_vector, n_aa_vector))
                                            if np.degrees(n_angle_hbond) >90.0 and np.degrees(n_angle_hbond)<180:
                                                rejection_list.append(n_atom)
                                        else:
                                            continue
# The internal hydrogen bond look up ends here. All the internally hydrogen bonded side chain dnor atoms are rejected.
                                if (n_atom not in rejection_list) and (n_atom != the_NE1_atom) and (atom_dic in donor_dict):
# Change this only if you need to. One representative atom from another side chain donor
#                                    CX = str("CB")
# Currently the CX has been set to CB atom. it could also be the iplus2 atom. Uncomment as per your need.
#                                        for na in n_res_atoms:
#                                            if (na.get_serial_number()) == iplus2:
#                                                CX = na.get_id()
#                                            else:
#                                                print("The Residue lacks i+2th Atom")
#                                    representative_atom = (n_atom.get_id())
# Pair of TRP - NE1 and Side Chain donor atoms are confirmed as the_NE1_atom and representative_atom.
# Carves the segment of query pdb id for CLICK alignment.
                                    carve(structure, model, chain, residue, n_atom, neighbour_atoms)
                                    for dataset_file in all_dataset_pdb:
                                        # Very important - format of pdb file is dataset/XXYY.pdb, hence int(16)
                                        if len(dataset_file) == int(16): #and (dataset_file == "dataset/1gbq.pdb"):
# Saves everything within click_output folder for CLICK alignment
                                            the_path = str(n_atom.get_full_id()[0])+"_"+str(n_atom.get_full_id()[1])+"_"+str(n_atom.get_full_id()[2])+"_"+str(the_NE1_atom.get_full_id()[3][1])+"_"+str(n_atom.get_full_id()[3][1])+"_"+str(n_atom.get_parent().get_resname()+"_"+n_atom.get_id()+"_"+str(dataset_file[-8:-4]))
# Mask Template PDB atoms for CLICK alignment
                                            mask_temp_Atoms(dataset_file, str(n_atom.get_id()), str(n_atom.get_parent().get_resname()), "_rnmd_ds", the_path)
# Mask Query Segment atoms which was "Carved" for CLICK alignment
                                            mask_query_Atoms (structure.get_id()+"__crvd__"+str(model.get_id())+chain.get_id()+str(residue.get_id()[1])+'_'+str(n_atom.get_parent().get_id()[1])+'_'+str(n_atom.get_id())+'.pdb', residue, n_atom.get_parent(), n_atom.get_id(), '__rnmd', the_path)


# Calls the above function.
#neighbour_search(input_structure)
