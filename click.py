#!/usr/bin/python
import sys
import os
import commands
import string

dna_atom=['  C',' DC',' DG','  G',' DA',' DT',' DI','  A','  U',' DU','  N','  I']	

alignThr = 1
saveThr  = 1
proteinA = sys.argv[1]
proteinB = sys.argv[2]

#Create proteinA and proteinB without DNA and HETATM atoms, and backup original proteinA and proteinB in proteinA.het and proteinB.het
os.system('mv '+proteinA+' '+proteinA+'.het')
fss = open(proteinA, 'w')	    
kt  = 0		     		
for line_atom in open(proteinA+'.het').readlines():
    if kt<2:
       if (line_atom[0:4]=='ATOM' and (line_atom[17:20] in dna_atom)==False):
          fss.write(line_atom)	
	  kt = 1
       if (line_atom[0:6]=='ENDMDL' and kt==1): break 			       	
fss.close()

os.system('mv '+proteinB+' '+proteinB+'.het')
fss = open(proteinB, 'w')	    
kt  = 0		     		
for line_atom in open(proteinB+'.het').readlines():
    if kt<2:
       if (line_atom[0:4]=='ATOM' and (line_atom[17:20] in dna_atom)==False):
          fss.write(line_atom)	
	  kt = 1
       if (line_atom[0:6]=='ENDMDL' and kt==1): break 			       	
fss.close()


#Use Modeller to get secondary structure and solvent accessibility values for proteinA and proteinB without DNA and HETATM atoms
command1='mod9.10 get_sec_struct.py '+proteinA+'  '+ proteinB
print command1
commands.getstatusoutput(command1)

#Write only secondary structure and solvent accessibility values into *.ss and *.sa files
command2='python get_struct.py '+proteinA+'  '+ proteinB
print command2
commands.getstatusoutput(command2)

#Move proteinA.het and proteinB.het back to proteinA and proteinB
os.system('mv '+proteinA+'.het '+proteinA)
os.system('mv '+proteinB+'.het '+proteinB)


#Run Click to find the best global alignment between two proteins based on the clique matching algorithm
command3='./click '+proteinA+'  '+ proteinB+' -a '+str(alignThr)+' -s '+str(saveThr)
print command3
commands.getstatusoutput(command3)

#Remove secondary structure and sovel accessibility files of proteinA and proteinB
if os.path.exists(proteinA+'.ss'):
   os.system('rm '+proteinA+'.s*')
   os.system('rm '+proteinA+'.psa')

if os.path.exists(proteinB+'.ss'):
   os.system('rm '+proteinB+'.s*')
   os.system('rm '+proteinB+'.psa')

