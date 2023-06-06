# This program analyzes a PDB structure to identify possible binding sites of the PPII helix along with the best PPII template.
# Input: 4-letter PDB ID (e.g., 1CKA)
# Output: If no binding site for PPII is found, a notification will be provided in the log file (pp2_pred_db_log.log).
#         If a binding site is found, the program will generate a "CLICK" transferred PPII on the binding site,
#         saved as input_pdb_best_template_model_pdb_4sim.pdb. Additionally, input_pdb_best_template_model_pdb_result.pdb
#         will contain all the accepted snapshots of the Monte Carlo moves.
# Author: Shashank Pritam (shashankpritam@gmail.com)
# License: LGPL

# Required Modules: numpy, scipy, modeller10.4 and biopython
# Python Version: 3.8.10
# Tested on: WSL Ubuntu 20.4
# Active internet connection is required if PDB files are not provided in the pp2pred folder (database_folder)

import sys
import os
import glob
import stat
import shutil
import logging
import numpy as np
import Bio.PDB
import warnings
from pathlib import Path
from Bio.PDB import PDBIO, PDBList
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

from modules.unmask_Atoms import unmask_Atoms_save
from modules.save_pdb import save_pdb
from modules.trp_nbr_lookup import neighbour_search


# Load the parameter file (param.txt), get the current working directory, and initialize the Biopython parser.
# Set up the log file.
parameter_file = open("param.txt").read().split('\n')
current_working_dir = os.getcwd()
parser = Bio.PDB.PDBParser(QUIET=True)
pdbl = PDBList()
io = PDBIO()
LOG_FILENAME = 'pp2_pred_db_log.log'
logging.basicConfig(filename=LOG_FILENAME, level=logging.INFO)


# Input, Output, and Command Handling
# If a database search is required, run pp2_pred_for_database.py and provide the query PDB ID in input_database.txt.
# Set the input command as sys.argv[1]. Alternatively, remove sys.argv[1] and uncomment the line below for interactive input.
input_pdb_given = sys.argv[1]  # input("Enter the four-letter PDB code of the query protein: ")
print("Input PDB Given: " + input_pdb_given)
# Provide the model of relevance here or define it as input.
input_receptor_model = 0

# Load the input PDB file using Biopython. If it is already present in the database_folder, it will be loaded instantly.
# Otherwise, Biopython will download the PDB from the RCSB website and load it into the module.
try:
    input_pdb = current_working_dir + "/database_folder/" + input_pdb_given + ".pdb"
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

# A function to get the peptide chain of a PDB from the template dataset containing 39 PDBs
def get_pep_chain(input_pdb):
    pdb_id_upper = input_pdb.upper()
    pdb_id_index = pdb_id.index(pdb_id_upper)
    input_peptide_chain = ppii_chain[pdb_id_index]
    return input_peptide_chain

# This information is needed to decide acceptor, donor, and acceptor antecedent atoms for hydrogen bond inference.
# Loads file with a list of acceptors, donors, and acceptor antecedents, creates a dictionary for acceptor_antecedent

# File with a list of acceptors
acceptor_file = open("H_bond_Acceptor_List").read().split('\n')
acceptor_file = filter(None, acceptor_file)
acceptor = []
acceptor_list = []
for i in acceptor_file:
    tmp = i.split('\t')
    acceptor.append(tmp[0] + ':' + tmp[1])
    acceptor_list.append(tmp[-1])

# File with a list of donors
donor_file = open("H_bond_Donor_List").read().split('\n')
donor_file = filter(None, donor_file)
donor = []
for i in donor_file:
    donor.append(i.replace('\t', ':'))

# Acceptor Antecedent Dictionary
acceptor_antecedent = {}
for atom in acceptor:
    acceptor_antecedent[atom] = acceptor_list[acceptor.index(atom)]

# A function to get the acceptor antecedent atom of an acceptor atom involved in a hydrogen bond
def aa_atom(atom):
    print(acceptor_list[acceptor.index(atom)])

# Subdirectories containing the template PDB files and the hydrogen bond info
all_dataset_pdb = glob.glob('dataset/*.pdb')
# This hydrogen bond info had 4 Angstrom as the distance cutoff between the acceptor and the donor atoms
# The range for the hydrogen bond angle was (90, 180)
hbond_files = glob.glob('data_hbond/*.txt')


# Look for a pair of TRP & NBR. Read neighhbor_look_up_4_pair_of_res for further information.
neighbour_search(input_structure)

# Set the CLICK folder.
files_4_click = glob.glob(f"{current_working_dir}/click_output/*/", recursive=True)

for folders in files_4_click:
    if folders.startswith(f"{current_working_dir}/click_output/{input_pdb_given.upper()}") or folders.startswith(f"{current_working_dir}/click_output/{input_pdb_given.lower()}"):
        renamed_pdb = glob.glob(f"{folders}/*_rnmd.pdb")
        dataset_renamed_file = glob.glob(f"{folders}/*_rnmd_ds.pdb")

        # Within all subfolders of this CLICK folder, a pair of query PDB and Template PDB is present for all (TRP, NBR_Atom) for all template PDBs
        # Call the function to structurally align these two PDBs
        # This function runs the CLICK alignment for each segment/PDB with all dataset PDB files
        # Before calling this function, make sure that the Parameters.inp file contains the desired atoms as representative_atoms.
        # Important to note that the representative_atoms in the first line of the CLICK parameter file Parameters.inp requires
        # an atom that occupies 4 character spaces including blank spaces. Please refer to http://cospi.iiserpune.ac.in/click/Contact/Contactus.jsp
        cmd = f'./click {renamed_pdb[0]} {dataset_renamed_file[0]} >/dev/null 2>&1'
        os.system(cmd)


# The atoms are masked with these "masks only"
click_atoms = ["AA", "BB", "CC", "NX", "EE"]
predicted_alignments = []

for folder in files_4_click:
    if folder.startswith(f"{current_working_dir}/click_output/{input_pdb_given.upper()}") or \
       folder.startswith(f"{current_working_dir}/click_output/{input_pdb_given.lower()}"):
        click_files = glob.glob(f"{folder}/*.clique")

        if len(click_files) > 0:
            with open(click_files[0], 'r') as infile:
                data = infile.read().split()
                matched_atoms = int(data[6])
                rmsd = float(data[9])
                structure_overlap = float(data[13])
                carved_frag_info = click_files[0].split('/')[-1]

                # Among all those aligned query and templated PDB with (CZ3, CA and NE1) of TRP and (CB and NX) of donor
                # Which ones are the "passable" alignments
                if matched_atoms == 5:  # RMSD <= 0.6 and
                    predicted_alignments.append(f"{carved_frag_info}_{rmsd}_{structure_overlap}")
                '''
                elif matched_atoms == 5:  # RMSD <= 1.0 and
                    predicted_alignments.append(f"{carved_frag_info}_{rmsd}_{structure_overlap}")
                '''
        else:
            print(f"No '.clique' files found in folder: {folder}")


# Remove self-alignments from passable alignments (e.g., 1CKA vs 1CKA)
unique_alignments = [
    alignment.split('_') for alignment in predicted_alignments 
    if alignment.split('_')[0].casefold() != alignment.split('_')[8][-4:].casefold()
]

filtered_alignments = []

# Replace the RMSD of CLICK with Total_RMSD = RMSD_TRP + RMSD_TRP_C + RMSD_TRP_O + RMSD_NBR_C
for alignment in unique_alignments:
    query_align = [alignment[0], alignment[4][0], alignment[4][1], alignment[4][2:], alignment[5], alignment[6], alignment[8][-4:]]

    for folder in files_4_click:
        if folder.startswith(f"{current_working_dir}/click_output/{input_pdb_given.upper()}") or folder.startswith(f"{current_working_dir}/click_output/{input_pdb_given.lower()}"):
            dataset_renamed_file = glob.glob(f"{folder}/*_rnmd.1.pdb")
            renamed_pdb = glob.glob(f"{folder}/*_rnmd_ds.1.pdb")
            click_file = glob.glob(f"{folder}/*.clique")

            fragment_info = click_file[0].split('/')[-1].split("_")
            fragment_info = [fragment_info[0], fragment_info[4][0], fragment_info[4][1], fragment_info[4][2:], fragment_info[5], fragment_info[6], fragment_info[8][5:9]]
            
            if fragment_info != query_align:
                continue

            old_predicted_receptor_structure = parser.get_structure(query_align[0], renamed_pdb[0])

            # Extract required atoms for alignment
            try:
                atom_data = {
                    "AA": old_predicted_receptor_structure[int(query_align[1])][query_align[2]][int(query_align[3])]["AA"],
                    "BB": old_predicted_receptor_structure[int(query_align[1])][query_align[2]][int(query_align[3])]["BB"],
                    "CC": old_predicted_receptor_structure[int(query_align[1])][query_align[2]][int(query_align[3])]["CC"],
                    "NX": old_predicted_receptor_structure[int(query_align[1])][query_align[2]][int(query_align[4])]["NX"],
                    "EE": old_predicted_receptor_structure[int(query_align[1])][query_align[2]][int(query_align[4])]["EE"],
                }
            except:
                continue

            old_template_peptide_structure = parser.get_structure(alignment[8][-4:], dataset_renamed_file[0])
            trp_trp_list = []
            nbr_nbr_list = []

            # Find the nearest TRP and NBR atoms
            for t_model in old_template_peptide_structure:
                for t_chain in t_model:
                    if str(t_chain.get_id()) == get_pep_chain(alignment[8][-4:]):
                        peptide_chain = t_chain
                    for t_residue in t_chain:
                        for t_atom in t_residue:
                            if t_atom.get_id() == "AA":
                                trp_trp_list.append([atom_data["AA"] - t_atom, t_atom])
                            elif t_atom.get_id() == "NX":
                                nbr_nbr_list.append([atom_data["NX"] - t_atom, t_atom])

            nearest_trp = sorted(trp_trp_list, key=lambda x: float(x[0]))[0]
            nearest_nbr = sorted(nbr_nbr_list, key=lambda x: float(x[0]))[0]

            query_trp_atoms = Bio.PDB.Selection.unfold_entities(atom_data["AA"].get_parent(), 'A')
            temp_trp_atoms = Bio.PDB.Selection.unfold_entities(nearest_trp[1].get_parent(), 'A')
            query_dnr_atoms = Bio.PDB.Selection.unfold_entities(atom_data["NX"].get_parent(), 'A')
            temp_dnr_atoms = Bio.PDB.Selection.unfold_entities(nearest_nbr[1].get_parent(), 'A')
            peptide_residues = Bio.PDB.Selection.unfold_entities(peptide_chain, 'R')
            n_num = 0

            # Count the number of neighboring atoms around the peptide residues
            for residue in peptide_residues:
                receptor_chain_atoms = Bio.PDB.Selection.unfold_entities(old_predicted_receptor_structure, 'A')
                neighbourhood_search_alpha_c = Bio.PDB.NeighborSearch(receptor_chain_atoms)
                try:
                    neighbour_atoms_alpha_c = neighbourhood_search_alpha_c.search(residue["CA"].coord, 3.5)
                except:
                    continue
                n_num += len(neighbour_atoms_alpha_c)

            trp = []
            nbr = []

            # Calculate RMSD for TRP and NBR atoms
            for q_atom, t_atom in zip(query_trp_atoms, temp_trp_atoms):
                if q_atom.get_id() == t_atom.get_id():
                    dist = q_atom - t_atom
                    trp.append(dist**2)

            for q_atom, t_atom in zip(query_dnr_atoms, temp_dnr_atoms):
                if q_atom.get_id() == t_atom.get_id():
                    rmsd_nbr = q_atom - t_atom
                    nbr.append(rmsd_nbr)

            RMSD_TRP = np.sqrt(np.sum(trp))
            RMSD_NBR = np.sqrt(np.sum(nbr))
            Total_RMSD = RMSD_TRP + RMSD_NBR
            Total_RMSD = f"{Total_RMSD:.2f}"
            align_details = (alignment[0], alignment[4][0], alignment[4][1], alignment[4][2:], alignment[5], alignment[6], alignment[8][-4:])
            score = f"{float(Total_RMSD) + float(n_num / 150):.2f}"
            new_table_data = [query_align, Total_RMSD, f"{RMSD_TRP:.2f}", f"{RMSD_NBR:.2f}", score, n_num]
            filtered_alignments.append(new_table_data)
            alignment[-2] = str(score)


# Sort the alignments based on score and RMSD(TRP)
for item in sorted(filtered_alignments, key=lambda x: (float(x[-2]), float(x[-4]))):
    print(item)

# If no "unique" alignments are found, report and exit. Otherwise, continue with the top alignment as the best alignment considered.
if filtered_alignments:
    best_alignment = sorted(filtered_alignments, key=lambda x: (float(x[-2]), float(x[-4])))[0]
    best_alignment_query = [best_alignment[0]]
else:
    logging.info(f"No binding site found in the given query structure {input_pdb}")
    sys.exit("No binding site found in the given query structure")

# If an alignment has Total_RMSD >= float(2.5), reject that alignment
if float(best_alignment[-2]) >= float(25000):
    logging.info(f"No binding site found in the given query structure {input_pdb}")
    sys.exit("No binding site found in the given query structure")

print(f"For query structure {best_alignment[0][0]}, predicted binding site details are - Model = {best_alignment[0][1]}, Chain = {best_alignment[0][2]}, TRP = {best_alignment[0][3]}, NBR = {best_alignment[0][4]}")
print(f"Template PPII is {best_alignment[0][6:]}, with a Score of {best_alignment[-2]}")

# Reports alignment with least RMSD with "passable" overlap.
#print(f"Best alignment details are - PDB ID = {best_alignment[8][-4:]}, RMSD = {best_alignment[-2]}, SO = {best_alignment[-1]}")
#print(f"{best_alignment[0]} {best_alignment[4][0]} {best_alignment[4][1]} {best_alignment[4][2:]}, {best_alignment[5]}, {best_alignment[8][-4:]}, {best_alignment[-2]}, {best_alignment[-1]}")

#print(best_alignment)

'''
global predicted_receptor_structure
global template_peptide_structure

# The below function returns the predicted receptor chain from the query structure and best aligned PPII chain from the template PDBs.
for folders in files_4_click:
    if folders.startswith(f"{current_working_dir}/click_output/{input_pdb_given.upper()}") or folders.startswith(f"{current_working_dir}/click_output/{input_pdb_given.lower()}"):
        dataset_renamed_file = glob.glob(f'{folders}/*_rnmd.1.pdb')
        renamed_pdb = glob.glob(f'{folders}/*_rnmd_ds.1.pdb')
        click_file = glob.glob(f'{folders}/*.clique')
        carved_frag_info = click_file[0].split('/')[-1].split("_")
        carved_frag_info = [carved_frag_info[0], carved_frag_info[4][0], carved_frag_info[4][1], carved_frag_info[4][2:], carved_frag_info[5], carved_frag_info[6], carved_frag_info[8][5:9]]
        
        if carved_frag_info == best_alignment_query[0]:
            receptor_chain = best_alignment[0][2]
            peptide_chain = get_pep_chain(best_alignment[0][6])
            unmask_Atoms_save(renamed_pdb[0], receptor_chain)
            unmask_Atoms_save(dataset_renamed_file[0], peptide_chain)
            
            renamed_pdb_unmask = f"{renamed_pdb[0][:-4]}_new.pdb"
            dataset_renamed_file_unmask = f"{dataset_renamed_file[0][:-4]}_new.pdb"
            
            save_pdb(renamed_pdb_unmask, best_alignment[0][0], receptor_chain)
            save_pdb(dataset_renamed_file_unmask, best_alignment[0][6], peptide_chain)
            
            predicted_receptor_structure = parser.get_structure(best_alignment[0][0], renamed_pdb_unmask)
            template_peptide_structure = parser.get_structure(best_alignment[0][6], dataset_renamed_file_unmask)
    #return (predicted_receptor_structure, receptor_chain, template_peptide_structure, peptide_chain)


#run_sim_sim (input_pdb_given)


input_receptor_chain_given = predicted_receptor_structure.get_id()
input_peptide_chain_given = template_peptide_structure.get_id()

input_receptor_chain = receptor_chain
input_peptide_chain = peptide_chain



input_receptor_structure = predicted_receptor_structure.get_id()
input_peptide_structure = template_peptide_structure.get_id()

# Below routines basically tidy ups the files for Monte-Carlo simulations
# Add the predicted_receptor_structure and template_peptide_structure onto a single file.
input_receptor_chain_given = run_sim_sim(input_pdb_given)[0].get_id()
input_peptide_chain_given = run_sim_sim(input_pdb_given)[2].get_id()

input_receptor_chain = run_sim_sim(input_pdb_given)[1]
input_peptide_chain = run_sim_sim(input_pdb_given)[3]



input_receptor_structure = run_sim_sim(input_pdb_given)[0]
input_peptide_structure = run_sim_sim(input_pdb_given)[2]


print("Receptor PDB ID : "+input_receptor_chain_given)
print("Peptide PDB ID : "+input_peptide_chain_given)
print(input_receptor_chain_given, input_receptor_chain, input_peptide_chain_given, input_peptide_chain)





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

result_filename = current_working_dir+"/"+str(input_receptor_chain_given)+str("_sim_result.pdb")
with open(result_filename, 'w+') as outfile:
    for files in glob.glob('pdb_traj_saved?.pdb'):
        with open(files) as file:
                for line in file:
                    outfile.write(line)
command12 = "sed -i 's/END/ENDMDL/g' "+ result_filename
#os.system(command12)


'''
#Important Snippet
#Remove All Residue Files - The Click Folder and MCMD PDBs
all_dump = [i for i in glob.glob("*.pdb") if "__crvd__" in i or "pdb_traj" in i]
for i in all_dump:
    os.remove(i)
    

def remove_read_only_files(func, path, excinfo):
    os.chmod(path, stat.S_IWRITE)
    func(path)

shutil.rmtree(current_working_dir+"/click_output", onerror=remove_read_only_files)
Path(current_working_dir+"/click_output").mkdir(parents=True, exist_ok=True)


#---------------------------------------------------------------End of The Line --------------------------------------------------------
