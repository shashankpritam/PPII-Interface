# Author - @Shashank Pritam - (shashankpritam@gmail.com).
# License - LGPL
# Required Modules -- biopython
# Working Python Version --  3.8.10 and tested system -- WSL Ubuntu 20.4
import os
import Bio.PDB
import glob
from Bio.PDB import PDBIO, Select, PDBList
current_working_dir = os.getcwd()
parser = Bio.PDB.PDBParser(QUIET=True)
pdbl = PDBList()
io = PDBIO()
test_pdb = glob.glob('test_pdb/*.pdb')
Main_Chain_Atoms = ["C", "N", "O", "CA"]
# The below function takes structure, model, chain, residue, n_atom, neighbour_atoms as input and carves that segment
# which is in the neighbourhood within some set cut-off (read param.txt) of NE1 of TRP
# To carve pdb make else condition return a value of 0. If it is 1 then it selects every thing.
# Simply put - You are getting a segment of a chain from the given structure which contains TRP and its neighbours within some
# set cut-off -neighbourhood_look_up_cut_off.
# Now each segment will have at least one another donor atom. These atoms are used for Click to find binding site.
# If Commented out; only predicted chain is taken considered for Click and no fragmentation of the PDB chain is taking place.
class ReceptorSelect(Select):
    def accept_model(self, model):
        if model.get_id() == int(input_receptor_model):
            return 1
        else:
            return 0
    def accept_residue(self, residue):
        if residue.id[0] == " ":
            return 1
        else:
            return 0

class PeptideSelect(Select):
    def accept_residue(self, residue):
        if residue.id[0] == " ":
            return 1
        else:
            return 0

class MCAtomSelect(Select):
    def accept_atom(self, atom):
        if atom.get_id() in Main_Chain_Atoms:
            return 1
        else:
            return 0
class AtomSelect(Select):
    def accept_atom(self, atom):
        if atom in neighbour_atoms:
            return 1
        else:
            return 0
#Class which selects only provided chain and none other
class ChainSelect(Select):
    def accept_chain(self, chain):
        if str(chain.get_id()) == str(input_chain.get_id()):
            return 1
        else:
            return 0

def carve(structure, input_model, input_chain, input_residue, n_atom, neighbour_atoms):
    pdb_id = structure.get_id()
    io.set_structure(structure)
# Remember to put AtomSelect() or ChainSelect() in io.save, according to your need.
    io.save(pdb_id+"__crvd__"+str(input_model.get_id())+input_chain.get_id()+str(input_residue.get_id()[1])+'_'+str(n_atom.get_parent().get_id()[1])+'_'+str(n_atom.get_id())+'.pdb')

# Test case - Uncomment below line
#print(test_pdb)
