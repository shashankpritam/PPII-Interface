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
