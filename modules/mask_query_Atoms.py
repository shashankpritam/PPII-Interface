import os
import Bio.PDB
from pathlib import Path
from Bio.PDB import PDBIO, Select, PDBList
parser = Bio.PDB.PDBParser(QUIET=True)
pdbl = PDBList()
io = PDBIO()
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
    for model in query_structure:
        for chain in model:
            for residue in chain:
                if residue == the_trp:
                    the_trp_residue = residue
                elif residue == the_nbr:
                    the_nbr_residue = residue
                    the_dn = the_nbr_residue[the_nbr_dnr]
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
