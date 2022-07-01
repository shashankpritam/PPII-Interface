# This program takes a pdb structure as input and return possible binding site of the PPII Helix alog with best PPII template.
# Input - 4 letter PDB ID. eg - 1CKA
# Output - If it has no binding site for PPII; You'll be notified; please check the log file - pp2_pred_db_log.log also.
#          If it has - You'll get to know the binding site and you'll get "CLICK" transferred PPII on the binding site
#          with input_pdb_best_template_model_pdb_4sim.pdb as file.
#          And input_pdb_best_template_model_pdb_result.pdb which will contain all the accepted snapshots of the Monte Carlo moves.
# Author - @Shashank Pritam - (shashankpritam@gmail.com). This experiment/script is an extension of earlier work done by
#          @Dr Nguyen Thanh Binh.
# License - LGPL

# Required Modules -- numpy, scipy, and biopython
# Working Python Version --  3.8.10 and tested system -- WSL Ubuntu 20.4
# Active internet connection is also required if PDB files are not provided in the pp2pred folder - database_folder
# All the required modules are imported here
import sys
import os
import glob
import shutil
import stat
import logging
from pathlib import Path
import numpy as np
import Bio.PDB
import warnings
from numpy.linalg import norm
from scipy.spatial.transform import Rotation as R
from Bio.PDB.vectors import calc_angle
from Bio.PDB import PDBIO, Select, PDBList
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

# The parameter file param.txt file loads here, we get the current working directory, and the biopython parser is loaded.
# And log file sets-up
parameter_file = open("param.txt").read().split('\n')
current_working_dir = os.getcwd()
parser = Bio.PDB.PDBParser(QUIET=True)
pdbl = PDBList()
io = PDBIO()
LOG_FILENAME = 'pp2_pred_db_log.log'
logging.basicConfig(filename=LOG_FILENAME, level=logging.INFO)


# Input Output and Command. Uncomment as per your requirements
# If database search is required please run pp2_pred_for_database.py and provide the query PDB ID in input_database.txt
# Let the input command be as sys.argv[1]. Also please provide the PDB files in database_folder for database look-up.
# In case of  running the script with each time query required - remove sys.argv[1] and uncomment the below line.
input_pdb_given = sys.argv[1]  #input("Enter the four letter PDB code of Query Protein: ")
pdbl = PDBList()
print("Input PDB Given : "+input_pdb_given)
#Provide Model of relevance here or define as input
input_receptor_model = 0

# Load the input pdb file onto Biopython. If already present in the database_folder - will be loaded instantly.
# Else biopython will download the PDB from the RCSB website and then load it onto the module.
try:
    input_pdb = current_working_dir+"/database_folder/"+input_pdb_given+".pdb"
    input_structure = parser.get_structure(input_pdb_given, input_pdb)
except:
    input_pdb = pdbl.retrieve_pdb_file(input_pdb_given, obsolete=False, pdir="database_folder", file_format="pdb", overwrite=False)
    input_structure = parser.get_structure(input_pdb_given, input_pdb)




# List of 39 Dataset PDB IDs - Templates and their receptor chain and peptide (ppii) chain ID.
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

# A function to get peptide chain of a PDB from the template dataset containing 39 PDBs
def get_pep_chain(input_pdb):
    pdb_id_upper = (input_pdb.upper())
    pdb_id_index = pdb_id.index(pdb_id_upper)
    input_peptide_chain = ppii_chain[pdb_id_index]
    return input_peptide_chain
# A function to get receptor chain of a PDB from the template dataset containing 39 PDBs
def get_rec_chain(input_pdb):
    pdb_id_upper = (input_pdb.upper())
    pdb_id_index = pdb_id.index(pdb_id_upper)
    input_receptor_chain = receptor_chain[pdb_id_index]
    return input_receptor_chain


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


# donor_dict is a dictionary which contains all the side donor atoms from their respective residues in a dictionary format.
# The information of donor_dict is taken from the Modeller software.
donor_dict = [('ARG', 'NE'), ('ARG', 'NH1'), ('ARG', 'NH2'), ('ASN', 'ND2'), ('ASX', 'ND2'), ('CYS', 'SG'), ('GLN', 'NE2'), ('GLX', 'NE2'), ('HIS', 'ND1'), ('HSE', 'NE2'), ('HSP', 'ND1'), ('HSP', 'NE2'), ('LYS', 'NZ'), ('SER', 'OG'), ('THR', 'OG1'), ('TRP', 'NE1'), ('TYR', 'OH')]


# Subdirectories containing the template PDB files and the hydrogen bond info.
all_dataset_pdb = glob.glob('dataset/*.pdb')
# This hydrogen bond info had 4 Angstrom as distance cut-off between the acceptor and the donor atoms.
# The range for hydrogen bond angle was (90, 180).
hbond_files = glob.glob('data_hbond/*.txt')
the_hbond_file = glob.glob('data_hbond/*hbond_trp_all.txt')

# This function runs CLICK alignment for each segment/pdb with all dataset pdb files
# Before calling this function makes sure that the Parameters.inp files contains the desired atoms as representative_atoms.
# Important to note down that that representative_atoms in the first line of the CLICK parameter file Parameters.inp requires
# an atom which acquires 4 character space including blank spaces. Please see http://cospi.iiserpune.ac.in/click/Contact/Contactus.jsp
def click4all(input_pdb1, input_pdb2):
    cmd = './click '+str(input_pdb1[0])+' '+str(input_pdb2[0])+''+'>/dev/null 2>&1'
    os.system(cmd)

# This function takes and input PDB ID from the template dataset and returns
# the list of TRP involved in Hydrogen Bond with the PPII - Only for the 39 PDB in dataset
# See comment above hbond_files
'''def hbond_trp(input_pdb):
    for file in hbond_files:
        if file.split('/')[1][0:4] == input_pdb:
            list_of_trp = []
            with open(file, "r") as infile:
                lines = infile.readlines()
                for line in lines:
                    data = line.split(' ')
                    if "TRP" in line:
                        trp_res_id = data[data.index("TRP")+1]
                        if trp_res_id not in list_of_trp:
                            list_of_trp.append(trp_res_id)
    return list_of_trp'''

def hbond_trp(input_pdb):
    with open("data_hbond/hbond_trp_all.txt", "r") as infile:
        lines = infile.readlines()
        list_of_trp = []
        for line in lines:
            #print(line)
            data = line.split(' ')
            #print(data[3])
            if input_pdb in line:
                print(input_pdb, line)
                if str("TRP") in data:
                    #print(data[1], line)
                    trp_res_id = data[data.index("TRP")+1]
                    if trp_res_id not in list_of_trp:
                        list_of_trp.append(trp_res_id)
    return list_of_trp




# This function takes input of template pdb id and an side chain donor atom (scda), side chain donor residue (scdr) of the scda
# and the suffix which acts as identifier as template file. The last parameter save_path is path where the ouput file is
# supposed to be saved.
# This function replaces the scda of scdr with NX and CB atom of scdr as EE.
# This function also replaces the NE1, CA and CZ3 of the Tryptophan residues (TRPs) from the trp_list acquired from hbond_trp
# as AA, BB and CC respectively.
# After replacing the atom names this functions saves the new PDB as per provided inputs.
# Suffix is used for renaming the output file.
# SCDA = Side Chain Donor ATOM
# SCDR = Side Chain Donor RESIDUE
def mask_temp_Atoms (input_pdb, scda, scdr, suffix, save_path):
    with open(input_pdb, "r") as infile:
        save_path = "click_output/"+save_path
        Path(save_path).mkdir(parents=True, exist_ok=True)
#        Use mkdir in case your system doesn't support save_path
#        os.mkdir(save_path)
        file_name = input_pdb[-8:-4]+str(suffix)+'.pdb'
        completeName = os.path.join(save_path, file_name)
        the_temp_chain = get_rec_chain(input_pdb[-8:-4])
        trp_list = (hbond_trp(input_pdb[-8:-4]))

        with open(completeName, 'w+') as outfile:
            lines = infile.readlines()
            for line in lines:
                if line.startswith("ATOM"):
                    atm_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain_id = line[21].strip()
                    res_seq = line[22:26].strip()
                    #if res_seq in trp_list:
                    if atm_name == "NE1" and res_name == "TRP":# and str(the_temp_chain) == str(chain_id):
                        line = line.replace(atm_name, "AA ", 1)
                        outfile.write(line)
                    elif atm_name == "CA" and res_name == "TRP":#  and str(the_temp_chain) == str(chain_id):
                        line = line.replace(atm_name, "BB", 1)
                        outfile.write(line)
                    elif atm_name == "CZ3" and res_name == "TRP":#  and str(the_temp_chain) == str(chain_id):
                        line = line.replace(atm_name, "CC ", 1)
                        outfile.write(line)
                        #else:
                            #outfile.write(line)
                    elif atm_name == scda and res_name == scdr:# and str(the_temp_chain) == str(chain_id):
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
                    elif atm_name == "CB" and res_name == scdr:# and str(the_temp_chain) == str(chain_id):
                        line = line.replace(atm_name, "EE", 1)
                        outfile.write(line)
                    else:
                        outfile.write(line)


# This function takes input of query pdb id and an Tryptophan residue as Biopython Residue Object - the_trp
# Other input parameters are - the_nbr; which the scdr but as biopython RESIDUE object and the_nbr_dnr -
# which the scda but as biopython ATOM object
# The suffix which acts as identifier as query file. The last parameter save_path is path where the ouput file is
# supposed to be saved.
# This function replaces the the_nbr with NX and CB atom of the_nbr_dnr as EE.
# This function also replaces the NE1, CA and CZ3 of the Tryptophan residues (TRPs) from the_trp
# as AA, BB and CC respectively.
# After replacing the atom names this functions saves the new PDB as per provided inputs.
# Suffix is used for renaming the output file.
# SCDA = Side Chain Donor ATOM
# SCDR = Side Chain Donor RESIDUE
def mask_query_Atoms (input_pdb, the_trp, the_nbr, the_nbr_dnr, suffix, save_path):
    save_path = "click_output/"+save_path
    Path(save_path).mkdir(parents=True, exist_ok=True)
#        os.mkdir(save_path)
    file_name = input_pdb[:-4]+str(suffix)+'.pdb'
    completeName = os.path.join(save_path, file_name)
    query_structure = parser.get_structure(input_pdb[0:4], input_pdb)
    the_dn_str = str()
    the_nbr_residue_int = int()
    for model in query_structure:
        for chain in model:
            for residue in chain:
                if residue == the_trp:
                    the_trp_residue = residue
                elif residue == the_nbr:
                    the_nbr_residue = residue
                    the_nbr_residue_int = int(the_nbr_residue.get_id()[1])
                    the_dn = the_nbr_residue[the_nbr_dnr]
                    the_dn_str = str(the_dn.get_id())
                    #the_ip2 = the_nbr_residue["CX"]
    with open(input_pdb, "r") as infile:
        with open(completeName, 'w+') as outfile:
             for line in infile:
                 if line.startswith("ATOM"):
                     atm_name = line[12:16].strip()
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
                     elif atm_name == the_dn_str and int(res_seq) == int(the_nbr_residue.get_id()[1]):
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
                     elif atm_name == "CB" and int(res_seq) == the_nbr_residue_int:
                         line = line.replace(atm_name, "EE", 1)
                         outfile.write(line)
                     else:
                         outfile.write(line)
                 else:
                     outfile.write(line)

# This function "undoes" what masking function above does.
# That is after the CLICK alignments all the desired files will have "normal"/unmasked file n_atom_residue
# with _new as suffix for identifier.
def unmask_Atoms_save (input_pdb, input_chain):
    with open(input_pdb, "r") as infile:
        output_file = input_pdb[:-4]+"_new.pdb"
        file_path = str(Path(input_pdb).resolve())
        the_dn = (file_path.split('/')[-2].split('_')[-2])
        with open(output_file, 'w+') as outfile:
             for line in infile:
                 if line.startswith("ATOM"):
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

# This function takes a pdb filename with complete path, it is 4 letter ID and a Chain ID as as input.
# This function then returns a new file in the parent directory with _4sim as suffix identifier.
# That PDB file should contain the pdb id provided with only the chain which was provided as an input parameter.
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



# The below function takes structure, model, chain, residue, n_atom, neighbour_atoms as input and carves that segment
# which is in the neighbourhood within some set cut-off (read param.txt) of NE1 of TRP
# To carve pdb make else condition return a value of 0. If it is 1 then it selects every thing.
# Simply put - You are getting a segment of a chain from the given structure which contains TRP and its neighbours within some
# set cut-off -neighbourhood_look_up_cut_off.
# Now each segment will have at least one another donor atom. These atoms are used for Click to find binding site.
# If Commented out; only predicted chain is taken considered for Click and no fragmentation of the PDB chain is taking place.
def carve(structure, input_model, input_chain, input_residue, n_atom, neighbour_atoms):
#    class AtomSelect(Select):
#        def accept_atom(self, atom):
#            if atom in neighbour_atoms:
#                return 1
#            else:
#                return 1
# Class which selects only provided chain and none other
    class ChainSelect(Select):
        def accept_chain(self, chain):
            if str(chain.get_id()) == str(input_chain.get_id()):
                return 1
            else:
                return 0
    pdb_id = structure.get_id()
    io.set_structure(structure)
# Remember to put AtomSelect() or ChainSelect() in io.save, according to your need.
    io.save(pdb_id+"__crvd__"+str(input_model.get_id())+input_chain.get_id()+str(input_residue.get_id()[1])+'_'+str(n_atom.get_parent().get_id()[1])+'_'+str(n_atom.get_id())+'.pdb')



# This functions below - neighbour_search looks for a side chain donor atom for hydrogen bond which are within 12 Angstrom from any (TRP, NE1)
# and considers only those which are not inter-hydrogen bound within the same chain.
# Please do not set/change the parameters here.
# Input is simply a biopython structure Object
# Neighbourhood Search Parameters
neighbourhood_look_up_cut_off = float(parameter_file[5].split( )[2])
H_Bond_Cut_Off = float(parameter_file[6].split( )[2])
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
                        chain_atoms  = Bio.PDB.Selection.unfold_entities(model, 'A')
                        neighbourhood_search = Bio.PDB.NeighborSearch(chain_atoms)
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
neighbour_search(input_structure)
# CLICK folder is set.
files_4_click = glob.glob(current_working_dir+"/click_output/*/", recursive = True)
for folders in files_4_click:
    if folders.startswith(current_working_dir+"/click_output/"+input_pdb_given.upper()) or folders.startswith(current_working_dir+"/click_output/"+input_pdb_given.lower()):
        renamed_pdb = glob.glob(folders+'/*_rnmd.pdb')
        dataset_renamed_file = glob.glob(folders+'/*_rnmd_ds.pdb')
# Within all subfolder of this CLICK folder a pair of query PDB and Template PDB is present for all (TRP, NBR_Atom) for all template PDBs
# Calls the function to structurally align these two PDBs
        click4all(renamed_pdb, dataset_renamed_file)

# The atoms are masked with these "masks only"
click_atoms = ["AA", "BB", "CC", "NX", "EE"]
predicted_alignments = []
for folders in files_4_click:
    if folders.startswith(current_working_dir+"/click_output/"+input_pdb_given.upper()) or folders.startswith(current_working_dir+"/click_output/"+input_pdb_given.lower()):
        click_file = glob.glob(folders+'/*.clique')
        with open(click_file[0], 'r') as infile:
            data_to_read = infile.read()
            data = data_to_read.split()
            split_data = data_to_read.split()
            Matched_Atoms = int(split_data[6])
            RMSD = float(split_data[9])
            Structure_Overlap = float(split_data[13])
            carved_frag_info = (click_file[0]).split('/')[-1]
# Among all those aligned query and templated PDB with (CZ3, CA and NE1) of TRP and (CB and NX) of donor
# Which ones are the "passable" alignments
            if Matched_Atoms == 5:#RMSD <= 0.6 and
                predicted_alignments.append(carved_frag_info+"_"+str(RMSD)+"_"+str(Structure_Overlap))
'''            elif Matched_Atoms == 5:#RMSD <= 1.0 and
                predicted_alignments.append(carved_frag_info+"_"+str(RMSD)+"_"+str(Structure_Overlap))'''

# Among these passable alignments, remove the trivial self alignments. (eg - 1CKA vs 1CKA)
list_of_unique_alignment = []
for alignments in predicted_alignments:
    alignments = alignments.split('_')
    if (alignments[0]).casefold() != (alignments[8][-4:]).casefold():
        list_of_unique_alignment.append(alignments)




# Code below replace the RMSD of CLICK with Total_RMSD = RMSD_TRP+RMSD_TRP_C+RMSD_TRP_O+RMSD_NBR_C
for align in list_of_unique_alignment:
    align_query = [align[0], align[4][0], align[4][1], align[4][2:], align[5], align[6], align[8][-4:]]
    for folders in files_4_click:
        if folders.startswith(current_working_dir+"/click_output/"+input_pdb_given.upper()) or folders.startswith(current_working_dir+"/click_output/"+input_pdb_given.lower()):
            dataset_renamed_file = glob.glob(folders+'/*_rnmd.1.pdb')
            renamed_pdb = glob.glob(folders+'/*_rnmd_ds.1.pdb')
            click_file = glob.glob(folders+'/*.clique')
            carved_frag_info = ((click_file[0]).split('/')[-1]).split("_")
            #carved_frag_info.pop(5)
            carved_frag_info = [carved_frag_info[0], carved_frag_info[4][0], carved_frag_info[4][1], carved_frag_info[4][2:], carved_frag_info[5], carved_frag_info[6], carved_frag_info[8][5:9]]
            #print(carved_frag_info, best_alignment_query)
            if carved_frag_info == align_query:
                #print(dataset_renamed_file, renamed_pdb, click_file)
                old_predicted_receptor_structure = parser.get_structure(align[0], renamed_pdb[0])
                try:
                    aa_atom = old_predicted_receptor_structure[int(align_query[1])][str(align_query[2])][int(align_query[3])]["AA"]
                    bb_atom = old_predicted_receptor_structure[int(align_query[1])][str(align_query[2])][int(align_query[3])]["BB"]
                    cc_atom = old_predicted_receptor_structure[int(align_query[1])][str(align_query[2])][int(align_query[3])]["CC"]
                    nx_atom = old_predicted_receptor_structure[int(align_query[1])][str(align_query[2])][int(align_query[4])]["NX"]
                    ee_atom = old_predicted_receptor_structure[int(align_query[1])][str(align_query[2])][int(align_query[4])]["EE"]
                except:
                    continue
#                print(aa_atom, bb_atom, cc_atom, nx_atom, ee_atom)
                old_template_peptide_structure = parser.get_structure((align[8][-4:]), dataset_renamed_file[0])
                trp_trp_list = []
                nbr_nbr_list = []
                for t_model in old_template_peptide_structure:
                    for t_chain in t_model:
                        if str(t_chain.get_id()) == (get_pep_chain(align[8][-4:])):
                            peptide_chain = t_chain
                        for t_residue in t_chain:
                            for t_atom in t_residue:
                                if t_atom.get_id() == "AA":
                                    to_append = [aa_atom-t_atom, t_atom]
                                    trp_trp_list.append(to_append)
                                elif t_atom.get_id() == "NX":
                                    to_append = [nx_atom-t_atom, t_atom]
                                    nbr_nbr_list.append(to_append)
                                else:
                                    continue
                the_nearest_trp = sorted(trp_trp_list, key = lambda x: float(x[0]))[0]
                the_nearest_nbr = sorted(nbr_nbr_list, key = lambda x: float(x[0]))[0]
                query_trp_atoms  = Bio.PDB.Selection.unfold_entities(aa_atom.get_parent(), 'A')
                temp_trp_atoms = Bio.PDB.Selection.unfold_entities(the_nearest_trp[1].get_parent(), 'A')
                query_dnr_atoms = Bio.PDB.Selection.unfold_entities(nx_atom.get_parent(), 'A')
                temp_dnr_atoms = Bio.PDB.Selection.unfold_entities(the_nearest_nbr[1].get_parent(), 'A')
                peptide_residues = Bio.PDB.Selection.unfold_entities(peptide_chain, 'R')
                n_num = 0
                for residue in peptide_residues:
                    receptor_chain_atoms  = Bio.PDB.Selection.unfold_entities(old_predicted_receptor_structure, 'A')
                    neighbourhood_search_alpha_c = Bio.PDB.NeighborSearch(receptor_chain_atoms)
                    try:
                        neighbour_atoms_alpha_c = neighbourhood_search_alpha_c.search(residue["CA"].coord, 3.5)
                    except:
                        continue
                    n_num = n_num+len(neighbour_atoms_alpha_c)
                trp = []
                nbr = []
                for q_atoms in query_trp_atoms:
                    for t_atoms in temp_trp_atoms:
                        if q_atoms.get_id() == t_atoms.get_id():
                            dist = q_atoms-t_atoms
                            trp.append(dist**2)
                for q_atoms in query_dnr_atoms:
                    for t_atoms in temp_dnr_atoms:
                        if q_atoms.get_id() == t_atoms.get_id():
                            rmsd_nbr = q_atoms-t_atoms
                            nbr.append(rmsd_nbr)
                RMSD_TRP = (np.sum(trp))**0.5
                RMSD_NBR = (np.sum(nbr))**0.5
                Total_RMSD = RMSD_TRP+RMSD_NBR
                Total_RMSD = "{:.2f}".format(Total_RMSD)
                align_details = (align[0], align[4][0], align[4][1], align[4][2:], align[5], align[6], align[8][-4:])
                score = "{:.2f}".format(float(Total_RMSD)+float(n_num/50))
                print(Total_RMSD, "{:.2f}".format(RMSD_TRP), "{:.2f}".format(RMSD_NBR), score, n_num, align_details)
                align[-2] = str(score)


# If no "unique" alignments are found reports and exists. Else continues with top alignment
# as best alignment considered. The top alignment is selected based on RMSD and structural overlap (SO) acts as tie break up.
if len(list_of_unique_alignment) != 0:
    best_alignment = sorted(list_of_unique_alignment, key = lambda x: (float(x[-2]), -float(x[-1])))[0]
    best_alignment_query = [best_alignment[0], best_alignment[4][0], best_alignment[4][1], best_alignment[4][2:], best_alignment[5], best_alignment[6], best_alignment[8][-4:]]
else:
    logging.info("No binding site found in the given query structure "+input_pdb)
    sys.exit("No binding site found in the given query structure")

print(list_of_unique_alignment)
# If an alignmnet has Total_RMSD >= float(2.5), that alignmnet should be rejected
if float(best_alignment[-2]) >= float(25000):
    logging.info("No binding site found in the given query structure "+input_pdb)
    sys.exit("No binding site found in the given query structure")


print("For query structure", best_alignment[0], "predicted binding site details are - ", "Model =", best_alignment[4][0],  "Chain =", best_alignment[4][1], "TRP =", best_alignment[4][2:], "NBR =", best_alignment[5])
print("Template PPII is", best_alignment[8][-4:],  "with a Score of", best_alignment[-2])

# Reports alignment with least RMSD with "passable" overlap.
#print("Best alignment details are - ", "PDB ID = ", best_alignment[8][-4:],  "RMSD = ", best_alignment[-2], "SO = ", best_alignment[-1])
#print(best_alignment[0], best_alignment[4][0], best_alignment[4][1], best_alignment[4][2:], best_alignment[5], best_alignment[8][-4:], best_alignment[-2], best_alignment[-1])

#print(best_alignment)


# The below function returns the predicted receptor chain from the query structure and best aligned PPII chain from the template PDBs.
def run_sim_sim (input_pdb):
    for folders in files_4_click:
        if folders.startswith(current_working_dir+"/click_output/"+input_pdb_given.upper()) or folders.startswith(current_working_dir+"/click_output/"+input_pdb_given.lower()):
            dataset_renamed_file = glob.glob(folders+'/*_rnmd.1.pdb')
            renamed_pdb = glob.glob(folders+'/*_rnmd_ds.1.pdb')
            click_file = glob.glob(folders+'/*.clique')
            carved_frag_info = ((click_file[0]).split('/')[-1]).split("_")
            #carved_frag_info.pop(5)
            #print(folders)
            carved_frag_info = [carved_frag_info[0], carved_frag_info[4][0], carved_frag_info[4][1], carved_frag_info[4][2:], carved_frag_info[5], carved_frag_info[6], carved_frag_info[8][5:9]]
            #print(carved_frag_info, best_alignment_query)
            if carved_frag_info == best_alignment_query:
                receptor_chain = best_alignment[4][1]
                peptide_chain = get_pep_chain(best_alignment[8][-4:])
                unmask_Atoms_save(renamed_pdb[0], receptor_chain)
                unmask_Atoms_save(dataset_renamed_file[0], peptide_chain)
                renamed_pdb_unmask = renamed_pdb[0][:-4]+"_new.pdb"
                dataset_renamed_file_unmask = dataset_renamed_file[0][:-4]+"_new.pdb"
                save_pdb (renamed_pdb_unmask, best_alignment[0], receptor_chain)
                save_pdb (dataset_renamed_file_unmask, best_alignment[8][-4:], peptide_chain)
                predicted_receptor_structure = parser.get_structure(best_alignment[0], renamed_pdb_unmask)
                template_peptide_structure = parser.get_structure((best_alignment[8][-4:]), dataset_renamed_file_unmask)
    return (predicted_receptor_structure, receptor_chain, template_peptide_structure, peptide_chain)


#run_sim_sim (input_pdb_given)


# Below routines basically tidy ups the files for Monte-Carlo simulations
# Add the predicted_receptor_structure and template_peptide_structure onto a single file.
input_receptor_chain_given = run_sim_sim(input_pdb_given)[0].get_id()
input_peptide_chain_given = run_sim_sim(input_pdb_given)[2].get_id()

input_receptor_chain = run_sim_sim(input_pdb_given)[1]
input_peptide_chain = run_sim_sim(input_pdb_given)[3]



input_receptor_structure = run_sim_sim(input_pdb_given)[0]
input_peptide_structure = run_sim_sim(input_pdb_given)[2]


#print("Receptor PDB ID : "+input_receptor_chain_given)
#print("Peptide PDB ID : "+input_peptide_chain_given)
#print(input_receptor_chain_given, input_receptor_chain, input_peptide_chain_given, input_peptide_chain)


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


#io.set_structure(input_receptor_structure)
#io.save(input_receptor_chain_given+"_receptor.pdb")
#io.set_structure(input_peptide_structure)
#io.save(input_peptide_chain_given+"_peptide.pdb")

rec_structure = parser.get_structure(input_receptor_chain_given, input_receptor_chain_given+"_4sim.pdb")
pep_structure = parser.get_structure(input_peptide_chain_given, input_peptide_chain_given+"_4sim.pdb")

for model_rec in rec_structure:
    for chain_rec in model_rec:
        chain_rec.id = 'X'
        io.set_structure(rec_structure)
        io.save(input_receptor_chain_given+"__4merge.pdb")

for model_pep in pep_structure:
    for chain_pep in model_pep:
        chain_pep.id = 'Y'
        io.set_structure(pep_structure)
        io.save(input_peptide_chain_given+"__4merge.pdb")

new_rec_structure = parser.get_structure(input_receptor_chain_given, input_receptor_chain_given+"__4merge.pdb")
new_pep_structure = parser.get_structure(input_peptide_chain_given, input_peptide_chain_given+"__4merge.pdb")

rec_file = input_receptor_chain_given+"__4merge.pdb"
pep_file = input_peptide_chain_given+"__4merge.pdb"

combined_file = input_receptor_chain_given+"__"+input_peptide_chain_given+"__4sim.pdb"

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
command7 = input_receptor_chain_given+"_4sim.pdb"
command8 = input_peptide_chain_given+"_4sim.pdb"
# Remove the files
os.remove(command4)
os.remove(command5)
os.remove(command7)
os.remove(command8)


# PDB for simulation is set.
simulation_pdb = parser.get_structure(input_receptor_chain_given[0:4], input_receptor_chain_given+"__"+input_peptide_chain_given+"__4sim.pdb")


# Calculates the energy function called Number of Hydrogen Bond Score
# NHBS Hydrogen bond distance
nhbs_h_bond_distance = float(parameter_file[10].split( )[2])
Main_Chain_Atoms = ["C", "N", "O", "CA"]
def NHBS(structure):
    for model in structure:
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
        for p_atoms in peptide_atoms:
            for r_atoms in receptor_atoms:
                if (p_atoms)-(r_atoms) <= nhbs_h_bond_distance:
                    peptide_residue = ((p_atoms.get_parent().get_resname()))
                    receptor_residue = ((r_atoms.get_parent().get_resname()))
                    peptide_id = (p_atoms.get_name()+":"+peptide_residue)
                    receptor_id = (r_atoms.get_name()+":"+receptor_residue)
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
        return (hb**2)


# Calculates the energy function called Restrained Hydrogen Bond Score
# RHBS Constraints
lower_bound_dist = float(parameter_file[14].split( )[2])
upper_bound_dist = float(parameter_file[15].split( )[2])
lower_bound_ang = float(parameter_file[16].split( )[2])
upper_bound_ang = float(parameter_file[17].split( )[2])
scale_distance_factor = int(parameter_file[18].split( )[2])
def RHBS(structure):
    for model in structure:
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
                        dist_score = dist_score + scale_distance_factor*(dist-lower_bound_dist)*(dist-upper_bound_dist)
                    else:
                        dist_score = dist_score + 0
                    if np.degrees(r_angle_hbond)>lower_bound_ang and np.degrees(r_angle_hbond)<upper_bound_ang:
                        ang_score = ang_score + (np.degrees(r_angle_hbond)-lower_bound_ang)*(np.degrees(r_angle_hbond)-upper_bound_ang)
                        ang_score = np.radians(ang_score)
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
                    else:
                        dist_score = dist_score + 0
                    if np.degrees(p_angle_hbond)>lower_bound_ang and np.degrees(p_angle_hbond)<upper_bound_ang:
                        ang_score = ang_score + (np.degrees(p_angle_hbond)-lower_bound_ang)*(np.degrees(p_angle_hbond)-upper_bound_ang)
                        ang_score = np.radians(ang_score)
                    else:
                        ang_score = ang_score + 0
                else:
                    continue
        return (dist_score+(ang_score))


# Calculates the energy function called Clash Scores
# CS Clash Constraints
clash_dist = float(parameter_file[22].split( )[2])
MCMC_Scale = int(parameter_file[23].split( )[2])
MCSC_Scale = int(parameter_file[24].split( )[2])
SCSC_Scale = int(parameter_file[25].split( )[2])
def CS(structure):
    for model in structure:
        for chain in model:
            if chain.get_id() == "X":
                r_chain = chain
            elif chain.get_id() == "Y":
                p_chain = chain
            else:
                print("something is not right")
        receptor_atoms  = Bio.PDB.Selection.unfold_entities(r_chain, 'A')
        peptide_atoms = Bio.PDB.Selection.unfold_entities(p_chain, 'A')
        MCMC = 0
        MCSC = 0
        SCSC = 0
        for p_atoms in peptide_atoms:
            for r_atoms in receptor_atoms:
                dist = p_atoms-r_atoms
                if (dist < clash_dist) and (r_atoms.get_id() in Main_Chain_Atoms) and (p_atoms.get_id() in Main_Chain_Atoms):
                    MCMC = MCMC+1
                elif (dist < clash_dist) and (r_atoms.get_id() not in Main_Chain_Atoms) and (p_atoms.get_id() in Main_Chain_Atoms):
                    MCSC = MCSC+1
                elif (dist < clash_dist) and (r_atoms.get_id() in Main_Chain_Atoms) and (p_atoms.get_id() not in Main_Chain_Atoms):
                    MCSC = MCSC+1
                elif (dist < clash_dist) and (r_atoms.get_id() not in Main_Chain_Atoms) and (p_atoms.get_id() not in Main_Chain_Atoms):
                    SCSC = SCSC+1
                else:
                    continue
        clash_score = (MCMC_Scale*(MCMC))+(MCSC_Scale*(MCSC**0.5))+(SCSC_Scale*(SCSC**0.25))
        return (clash_score)


# Calculates the energy function called Pseudo-Van-der-Waal energy
# PVWS Constraints
PVWS_upper = float(parameter_file[28].split( )[2])
PVWS_lower = clash_dist
def PVWS(structure):
    for model in structure:
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
        for p_atoms in peptide_atoms:
            for r_atoms in receptor_atoms:
                dist = p_atoms-r_atoms
                if (dist > clash_dist) and (dist < PVWS_upper):
                    PVWS_num = PVWS_num+1
                else:
                    PVWS_num = PVWS_num+0
        return(PVWS_num)


# Calculates the total energy function
def total_energy(simulation_pdb):
    energy = -(NHBS(simulation_pdb))+(RHBS(simulation_pdb))+(CS(simulation_pdb))-(PVWS(simulation_pdb))
    return(energy)
#print("Total Energy =", total_energy(simulation_pdb), ", NHBS = ", (NHBS(simulation_pdb)), ", RHBS = ", (RHBS(simulation_pdb)), ", CS = ", (CS(simulation_pdb)), ", PVWS = ", (PVWS(simulation_pdb)))


#Monte Carlo Parameters
temperature = float(parameter_file[33].split( )[2])
#float(673.15) In Kelvin
gas_constant = float(parameter_file[36].split( )[2])
#In J⋅K−1⋅mol−1
translation_scale = int(parameter_file[39].split( )[2])
#translation =  translation_scale*(±0.5)
rotation_degrees = float(parameter_file[42].split( )[2])
#± rotation_degrees°
#Number of accepted moves you want the Monte Carlo Iteration
MC_accepted_move_limit = int(parameter_file[46].split( )[2])
#Number of accepted moves you want the Monte Carlo Iteration
MC_rejected_move_limit = int(parameter_file[47].split( )[2])

# This function calculates the metropolis parameter P
def metropolis_p(E1,E2):
    diff = float(E1-E2)
    div = float(gas_constant*temperature)
    diff_by_div = diff/div
    P = np.exp(diff_by_div)
    return P

# These transform functions or trans_x/y/z function transform the co-ordinates of the provided atom
# atom is a biopython atom object and dx/dy/dz is the amount by which we want to move the atom in
# respective directions
# the first line of each such function only translates the coordinates
# while the second one also rotates the coordinates by multiplying the coordinate vector
# by the rotation matrix - rot
# We can comment and uncomment these two lines according to our requirements
def x_trans(atom, dx, rot):
    #return atom.coord+(np.array([float(dx), 0, 0], dtype=float))
    return np.dot(rot, atom.coord)+(np.array([float(dx), 0, 0], dtype=float))
def y_trans(atom, dy, rot):
    #return atom.coord+(np.array([0, float(dy), 0], dtype=float))
    return np.dot(atom.coord, rot)+(np.array([0, float(dy), 0], dtype=float))
def z_trans(atom, dz, rot):
    #return atom.coord+(np.array([0, 0, float(dz)], dtype=float))
    return np.dot(atom.coord, rot)+(np.array([0, 0, float(dz)], dtype=float))



# This function move_pdb moves Y chain the provided structure
def move_pdb(structure, iteration):
    for model in structure:
        for chain in model:
            if chain.get_id() == "Y":
# For only the Y chain of the provided structure
# the_end_list contains the atom coordinates of atoms at the end of the chain Y
                the_end_list = np.empty((0,3), float)
                the_peptide_atoms = Bio.PDB.Selection.unfold_entities(chain, 'A')
                for atom in the_peptide_atoms:
                    atom_coord = np.array(atom.coord)
                    the_end_list = np.append(the_end_list, [atom_coord], axis=0)
                the_first_end = the_end_list[0]
                the_last_end = the_end_list[len(the_end_list)-1]
# mean_coordinate is the mean coordinate through all the peptide atoms
                mean_coordinate = (np.sum(the_end_list, axis=0)/len(the_end_list))
# direction_vector is the direction vector generated by the atoms at the both ends of the Y chain
                direction_vector = (the_first_end-the_last_end)
                direction_vector = direction_vector/ np.sqrt(np.sum(direction_vector**2))
                random_vector = mean_coordinate + ((np.random.uniform(-1, 1))*direction_vector)
                random_vector = random_vector/norm(random_vector)
# the random_vector is a random unit vector having the direction same as direction_vector
                rotn = R.from_rotvec(rotation_degrees*random_vector)
                rot = rotn.as_matrix()
# with the help of random_vector we generat a random matrix rot
                dx = translation_scale*(np.random.uniform(-0.5, 0.5))
                dy = translation_scale*(np.random.uniform(-0.5, 0.5))
                dz = translation_scale*(np.random.uniform(-0.5, 0.5))
                for residue in chain:
                    for atom in residue:
                        atom.coord = (x_trans(atom, dx, rot))
                        atom.coord = (y_trans(atom, dy, rot))
                        atom.coord = (z_trans(atom, dz, rot))
        io.set_structure(structure)
        file_name = "pdb_traj"+str(iteration)+".pdb"
        io.save(file_name)





energy_data = [0]
energy_plot = []
iteration_axis = []
# The function decide_move moves the coordinates of the peptide chain
# move_limit is a parameter which counts total moves but not the accepted moves; total rejected moves.
# if total rejected moves are higher than move_limit the Monte Carlo Simulations will stop.
def decide_move(input_structure, iteration, move_number):
# Limit on total accepted Monte Carlo moves = limit = MC_iteration
# Limit on total rejected Monte Carlo moves = move_limit
# E2_Val = total energy of current Monte Carlo Move
# While E1_Val is the total energy of the previous Monte Carlo Move
    limit = MC_accepted_move_limit
    move_limit = MC_rejected_move_limit
    E2_Val = float(total_energy(input_structure))
    last_energy_position = (len(energy_data)-1)
    E1_Val = float(energy_data[last_energy_position])
    while iteration <= limit and move_number <= move_limit:
# While the total accepted Monte Carlo move is less than provided MC_iteration
# and total rejected Monte Carlo move is less than provided move_limit
# Initialization of Monte Carlo
        if iteration == 0:
            energy_data[0] = (E2_Val)
            iteration = iteration+1
            move_pdb(input_structure, iteration)
        else:
            traj_file_name = "pdb_traj"+str(iteration)+".pdb"
            traj_file_str = parser.get_structure(input_receptor_chain_given[0:4], traj_file_name)
            saved_file_name = "pdb_traj_saved"+str(iteration)+".pdb"
            P_Val = (metropolis_p(E1_Val, E2_Val))
            if P_Val > 1:
# Accepted Move of Monte Carlo
                print("Accepted,", "p :", P_Val, "E :", E2_Val, "Total Accepted Moves :", iteration, "Total Rejected Moves :", move_number)
                energy_data.append(E2_Val)
                io.set_structure(input_structure)
                io.save(saved_file_name)
                energy_plot.append(E2_Val)
                iteration_axis.append(iteration)
                iteration = iteration+1
                move_pdb(traj_file_str, iteration)
                decide_move(traj_file_str, iteration, move_number)
            else:
# Rejected Move of Monte Carlo
                print("Rejected,", "p :", P_Val, "E :", E2_Val, "Total Accepted Moves :", iteration, "Total Rejected Moves :", move_number)
                iteration = iteration+0
                move_number = move_number + 1
                move_pdb(input_structure, iteration)
                decide_move(input_structure, iteration, move_number)
            break

# Monte Carlo Begins/Initializes here, comment the below line to just get the prediction site for ppii
# and avoid Monte Carlo
#decide_move(simulation_pdb, 0, 0)





##To create a single multi-model .pdb file - out_combined.pdb
result_filename = current_working_dir+"/"+str(input_receptor_chain_given)+str("_sim_result.pdb")
with open(result_filename, 'w+') as outfile:
    for files in glob.glob('pdb_traj_saved?.pdb'):
        with open(files) as file:
                for line in file:
                    outfile.write(line)
command12 = "sed -i 's/END/ENDMDL/g' "+ result_filename
os.system(command12)


#Important Snippet
#Remove All Residue Files - The Click Folder and MCMD PDBs
all_dump = [i for i in glob.glob("*.pdb") if "crvd" in i or "pdb_traj" in i]
for i in all_dump:
    os.remove(i)

def remove_read_only_files(func, path, excinfo):
    os.chmod(path, stat.S_IWRITE)
    func(path)

shutil.rmtree(current_working_dir+"/click_output", onerror=remove_read_only_files)
Path(current_working_dir+"/click_output").mkdir(parents=True, exist_ok=True)


#---------------------------------------------------------------End of The Line --------------------------------------------------------
