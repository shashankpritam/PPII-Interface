import glob
pdb = []
trp = []
c = []
o = []
ne1 = []
ca = []
nx = []
TRP_NBR = []

hbond_files = glob.glob('data_hbond/*.txt')

def hbond_trp(input_pdb, trp):
    for file in hbond_files:
        if file.split('/')[1][0:4] == input_pdb:
            list_of_trp = []
            infile = open(file, "r")
            lines = infile.read().split()
            if trp in lines:
                print("True", trp)
            else:
                print("False", trp)


with open ("trp_nbr_pred.txt", 'r') as infile:
    for line in infile:
        data = line.split(',')
        input_pdb = data[0]
        trp = data[1]
        hbond_trp(input_pdb, trp)











'''with open ("rmsd_result.txt", 'r') as infile:
    for line in infile:
        if line.startswith("Input PDB Given"):
            pdb.append(line.split(' ')[4])
            #print(line.split(' ')[4])
        elif line.startswith("For query structure"):
            TRP_NBR.append(line)
            print(line.split(' ')[4], line.split(' ')[-5], line.split(' ')[-1])
        elif line.startswith("RMSD_TRP"):
            trp.append(line.split(' ')[1])
            c.append(line.split(' ')[3])
            o.append(line.split(' ')[5])
            ne1.append(line.split(' ')[7])
            ca.append(line.split(' ')[9])
            nx.append(line.split(' ')[11])
            #print(line.split(' ')[0], line.split(' ')[1], line.split(' ')[2], line.split(' ')[3], line.split(' ')[4], line.split(' ')[5], line.split(' ')[6], line.split(' ')[7], line.split(' ')[8], line.split(' ')[9], line.split(' ')[10], line.split(' ')[11])
        else:
            continue
'''
'''with open ("Analysis of RMSD.txt", 'r') as infile:
    for line in infile:
        data = line.split(',')
        #print(data)
        if data[1] == data[7]: #== data[4]:#
            print("Yes")
        else:
            print("No")'''
