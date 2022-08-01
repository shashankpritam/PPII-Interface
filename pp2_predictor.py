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

from modules.hydrogen_bonds import hbond_trp
from modules.mask_temp_Atoms import mask_temp_Atoms
from modules.mask_query_Atoms import mask_query_Atoms
from modules.unmask_Atoms import unmask_Atoms_save
from modules.save_pdb import save_pdb
from modules.carve_pdb import carve
from modules.trp_nbr_lookup import neighbour_search
from modules.pep_simulation import decide_move



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
ppii_chain = []
with open("pdb_chains.txt", "r") as chain_info:
    next(chain_info)
    count = 0
    for line in chain_info:
        count += 1
        pdb_id.append(line.split()[0])
        ppii_chain.append(line.split()[2])

# A function to get peptide chain of a PDB from the template dataset containing 39 PDBs
def get_pep_chain(input_pdb):
    pdb_id_upper = (input_pdb.upper())
    pdb_id_index = pdb_id.index(pdb_id_upper)
    input_peptide_chain = ppii_chain[pdb_id_index]
    return input_peptide_chain




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




# Subdirectories containing the template PDB files and the hydrogen bond info.
all_dataset_pdb = glob.glob('dataset/*.pdb')
# This hydrogen bond info had 4 Angstrom as distance cut-off between the acceptor and the donor atoms.
# The range for hydrogen bond angle was (90, 180).
hbond_files = glob.glob('data_hbond/*.txt')


# This function runs CLICK alignment for each segment/pdb with all dataset pdb files
# Before calling this function makes sure that the Parameters.inp files contains the desired atoms as representative_atoms.
# Important to note down that that representative_atoms in the first line of the CLICK parameter file Parameters.inp requires
# an atom which acquires 4 character space including blank spaces. Please see http://cospi.iiserpune.ac.in/click/Contact/Contactus.jsp
def click4all(input_pdb1, input_pdb2):
    cmd = './click '+str(input_pdb1[0])+' '+str(input_pdb2[0])+''+'>/dev/null 2>&1'
    os.system(cmd)

# Look for pair of TRP & NBR. Read neighhbor_look_up_4_pair_of_res for further information.
#neighbour_search(input_structure)
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

#print(predicted_alignments)
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
# Monte Carlo Begins/Initializes here, comment the below line to just get the prediction site for ppii
# and avoid Monte Carlo
#decide_move(simulation_pdb, 0, 0)


##To create a single multi-model .pdb file - out_combined.pdb
'''result_filename = current_working_dir+"/"+str(input_receptor_chain_given)+str("_sim_result.pdb")
with open(result_filename, 'w+') as outfile:
    for files in glob.glob('pdb_traj_saved?.pdb'):
        with open(files) as file:
                for line in file:
                    outfile.write(line)
command12 = "sed -i 's/END/ENDMDL/g' "+ result_filename'''
#os.system(command12)


#Important Snippet
#Remove All Residue Files - The Click Folder and MCMD PDBs
all_dump = [i for i in glob.glob("*.pdb") if "__crvd__" in i or "pdb_traj" in i]
for i in all_dump:
    os.remove(i)

def remove_read_only_files(func, path, excinfo):
    os.chmod(path, stat.S_IWRITE)
    func(path)

#shutil.rmtree(current_working_dir+"/click_output", onerror=remove_read_only_files)
#Path(current_working_dir+"/click_output").mkdir(parents=True, exist_ok=True)


#---------------------------------------------------------------End of The Line --------------------------------------------------------
