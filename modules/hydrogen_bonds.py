# Author - @Shashank Pritam - (shashankpritam@gmail.com).
# License - LGPL
# Working Python Version --  3.8.10 and tested system -- WSL Ubuntu 20.4
# This function takes and input PDB ID from the template dataset and returns
# the list of TRP involved in Hydrogen Bond with the PPII - Only for the 39 PDB in dataset
# See comment above hbond_files
def hbond_trp(input_pdb):
    with open("data_hbond/hbond_trp_all.txt", "r") as infile:
        lines = infile.readlines()
        list_of_trp = []
        for line in lines:
            data = line.split(' ')
            if input_pdb.upper() in line:
                trp_res_id = str(data[8])
                chain_id = str(data[9][1])
                #print(trp_res_id, chain_id)
                res_chain_serial = [trp_res_id+chain_id]
                if res_chain_serial not in list_of_trp:
                    list_of_trp.append(res_chain_serial)
    return list_of_trp
