import os
from pathlib import Path
from modules.hydrogen_bonds import hbond_trp
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
#        os.mkdir(save_path, exist_ok=True)
        file_name = input_pdb[-8:-4]+str(suffix)+'.pdb'
        completeName = os.path.join(save_path, file_name)
        the_temp_trp_chain = (hbond_trp(input_pdb[-8:-4]))
        trp_list = (hbond_trp(input_pdb[-8:-4]))
        with open(completeName, 'w+') as outfile:
            lines = infile.readlines()
            for line in lines:
                if line.startswith("ATOM"):
                    atm_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain_id = line[21].strip()
                    res_seq = line[22:26].strip()
                    if [res_seq+","+chain_id] in hbond_trp(input_pdb[-8:-4]):
                        if atm_name == "NE1":
                            line = line.replace(atm_name, "AA ", 1)
                            outfile.write(line)
                        elif atm_name == "CA":
                            line = line.replace(atm_name, "BB", 1)
                            outfile.write(line)
                        elif atm_name == "CZ3":
                            line = line.replace(atm_name, "CC ", 1)
                            outfile.write(line)
                        else:
                            outfile.write(line)
                    elif atm_name == scda and res_name == scdr:
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
                    elif atm_name == "CB" and res_name == scdr:
                        line = line.replace(atm_name, "EE", 1)
                        outfile.write(line)
                    else:
                        outfile.write(line)
