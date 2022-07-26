import Bio.PDB
parameter_file = open("param.txt").read().split('\n')
Main_Chain_Atoms = ["C", "N", "O", "CA"]
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
