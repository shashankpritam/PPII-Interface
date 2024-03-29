# This program takes a pdb structure as input and return possible binding site of the PPII Helix alog with best PPII template.
# Author - @Shashank Pritam - (shashankpritam@gmail.com).
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

parser = Bio.PDB.PDBParser(QUIET=True)
pdbl = PDBList()
io = PDBIO()


# The parameter file param.txt file loads here, we get the current working directory, and the biopython parser is loaded.
# And log file sets-up
parameter_file = open("param.txt").read().split('\n')
from modules.energy_nhbs import NHBS
from modules.energy_rhbs import RHBS
from modules.energy_pvws import PVWS
from modules.energy_clash_score import CS

# Calculates the total energy function
def total_energy(simulation_pdb):
    energy = -(NHBS(simulation_pdb))+(RHBS(simulation_pdb))+(CS(simulation_pdb))-(PVWS(simulation_pdb))
    return(energy)
#print("Total Energy =", total_energy(simulation_pdb), ", NHBS = ", (NHBS(simulation_pdb)), ", RHBS = ", (RHBS(simulation_pdb)), ", CS = ", (CS(simulation_pdb)), ", PVWS = ", (PVWS(simulation_pdb)))


#Monte Carlo Parameters
temperature = float(parameter_file[33].split( )[2])
#float(673.15) In Kelvin
gas_constant = float(parameter_file[36].split( )[2])
#In J⋅K−1⋅mol−1
translation_scale = int(parameter_file[39].split( )[2])
#translation =  translation_scale*(±0.5)
rotation_degrees = float(parameter_file[42].split( )[2])
#± rotation_degrees°
#Number of accepted moves you want the Monte Carlo Iteration
MC_accepted_move_limit = int(parameter_file[46].split( )[2])
#Number of accepted moves you want the Monte Carlo Iteration
MC_rejected_move_limit = int(parameter_file[47].split( )[2])

# This function calculates the metropolis parameter P
def metropolis_p(E1,E2):
    diff = float(E1-E2)
    div = float(gas_constant*temperature)
    diff_by_div = diff/div
    P = np.exp(diff_by_div)
    return P

# These transform functions or trans_x/y/z function transform the co-ordinates of the provided atom
# atom is a biopython atom object and dx/dy/dz is the amount by which we want to move the atom in
# respective directions
# the first line of each such function only translates the coordinates
# while the second one also rotates the coordinates by multiplying the coordinate vector
# by the rotation matrix - rot
# We can comment and uncomment these two lines according to our requirements
def x_trans(atom, dx, rot):
    #return atom.coord+(np.array([float(dx), 0, 0], dtype=float))
    return np.dot(rot, atom.coord)+(np.array([float(dx), 0, 0], dtype=float))
def y_trans(atom, dy, rot):
    #return atom.coord+(np.array([0, float(dy), 0], dtype=float))
    return np.dot(atom.coord, rot)+(np.array([0, float(dy), 0], dtype=float))
def z_trans(atom, dz, rot):
    #return atom.coord+(np.array([0, 0, float(dz)], dtype=float))
    return np.dot(atom.coord, rot)+(np.array([0, 0, float(dz)], dtype=float))



# This function move_pdb moves Y chain the provided structure
def move_pdb(structure, iteration):
    for model in structure:
        for chain in model:
            if chain.get_id() == "Y":
# For only the Y chain of the provided structure
# the_end_list contains the atom coordinates of atoms at the end of the chain Y
                the_end_list = np.empty((0,3), float)
                the_peptide_atoms = Bio.PDB.Selection.unfold_entities(chain, 'A')
                for atom in the_peptide_atoms:
                    atom_coord = np.array(atom.coord)
                    the_end_list = np.append(the_end_list, [atom_coord], axis=0)
                the_first_end = the_end_list[0]
                the_last_end = the_end_list[len(the_end_list)-1]
# mean_coordinate is the mean coordinate through all the peptide atoms
                mean_coordinate = (np.sum(the_end_list, axis=0)/len(the_end_list))
# direction_vector is the direction vector generated by the atoms at the both ends of the Y chain
                direction_vector = (the_first_end-the_last_end)
                direction_vector = direction_vector/ np.sqrt(np.sum(direction_vector**2))
                random_vector = mean_coordinate + ((np.random.uniform(-1, 1))*direction_vector)
                random_vector = random_vector/norm(random_vector)
# the random_vector is a random unit vector having the direction same as direction_vector
                rotn = R.from_rotvec(rotation_degrees*random_vector)
                rot = rotn.as_matrix()
# with the help of random_vector we generat a random matrix rot
                dx = translation_scale*(np.random.uniform(-0.5, 0.5))
                dy = translation_scale*(np.random.uniform(-0.5, 0.5))
                dz = translation_scale*(np.random.uniform(-0.5, 0.5))
                for residue in chain:
                    for atom in residue:
                        atom.coord = (x_trans(atom, dx, rot))
                        atom.coord = (y_trans(atom, dy, rot))
                        atom.coord = (z_trans(atom, dz, rot))
        io.set_structure(structure)
        file_name = "pdb_traj"+str(iteration)+".pdb"
        io.save(file_name)





energy_data = [0]
energy_plot = []
iteration_axis = []
# The function decide_move moves the coordinates of the peptide chain
# move_limit is a parameter which counts total moves but not the accepted moves; total rejected moves.
# if total rejected moves are higher than move_limit the Monte Carlo Simulations will stop.
def decide_move(input_structure, iteration, move_number):
# Limit on total accepted Monte Carlo moves = limit = MC_iteration
# Limit on total rejected Monte Carlo moves = move_limit
# E2_Val = total energy of current Monte Carlo Move
# While E1_Val is the total energy of the previous Monte Carlo Move
    limit = MC_accepted_move_limit
    move_limit = MC_rejected_move_limit
    E2_Val = float(total_energy(input_structure))
    last_energy_position = (len(energy_data)-1)
    E1_Val = float(energy_data[last_energy_position])
    while iteration <= limit and move_number <= move_limit:
# While the total accepted Monte Carlo move is less than provided MC_iteration
# and total rejected Monte Carlo move is less than provided move_limit
# Initialization of Monte Carlo
        if iteration == 0:
            energy_data[0] = (E2_Val)
            iteration = iteration+1
            move_pdb(input_structure, iteration)
        else:
            traj_file_name = "pdb_traj"+str(iteration)+".pdb"
            traj_file_str = parser.get_structure(input_receptor_chain_given[0:4], traj_file_name)
            saved_file_name = "pdb_traj_saved"+str(iteration)+".pdb"
            P_Val = (metropolis_p(E1_Val, E2_Val))
            if P_Val > 1:
# Accepted Move of Monte Carlo
                print("Accepted,", "p :", P_Val, "E :", E2_Val, "Total Accepted Moves :", iteration, "Total Rejected Moves :", move_number)
                energy_data.append(E2_Val)
                io.set_structure(input_structure)
                io.save(saved_file_name)
                energy_plot.append(E2_Val)
                iteration_axis.append(iteration)
                iteration = iteration+1
                move_pdb(traj_file_str, iteration)
                decide_move(traj_file_str, iteration, move_number)
            else:
# Rejected Move of Monte Carlo
                print("Rejected,", "p :", P_Val, "E :", E2_Val, "Total Accepted Moves :", iteration, "Total Rejected Moves :", move_number)
                iteration = iteration+0
                move_number = move_number + 1
                move_pdb(input_structure, iteration)
                decide_move(input_structure, iteration, move_number)
            break

# Monte Carlo Begins/Initializes here, comment the below line to just get the prediction site for ppii
# and avoid Monte Carlo
