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
