from pathlib import Path
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
