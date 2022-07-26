import Bio.PDB
import numpy as np
from Bio.PDB.vectors import calc_angle
parameter_file = open("param.txt").read().split('\n')
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
