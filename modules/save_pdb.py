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
