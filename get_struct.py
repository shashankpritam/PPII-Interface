#!/usr/bin/python
import sys
import string
import os

#Write only secodary structure and sovel accessibility values from *.ssm and *.psa into *.ss and *.sa files
AA={'GLY': 'G','ALA': 'A','VAL': 'V','LEU': 'L','ILE': 'I','MET': 'M','PRO': 'P','PHE': 'F','TRP': 'W','TYR': 'Y','SER': 'S','THR': 'T','CYS': 'C','ASN': 'N','GLN': 'Q','LYS': 'K','ARG': 'R','HIS': 'H','ASP': 'D','GLU': 'E','XLE': 'J','GLX': 'Z','ASX': 'B','UNK': 'X','XAA': 'X'}

if os.path.exists(sys.argv[1]+'.ssm')==False: sys.exit(1)
if os.path.exists(sys.argv[1]+'.psa')==False: sys.exit(1)
if os.path.exists(sys.argv[2]+'.ssm')==False: sys.exit(1)
if os.path.exists(sys.argv[2]+'.psa')==False: sys.exit(1)

count = 0
fss=open(sys.argv[1]+'.ss', 'w')
for line in open(sys.argv[1]+'.ssm').readlines():
  count=count+1
  if count>13:
     cols = string.split(line)
     if (cols[2] in AA.keys()): fss.write(cols[-1]+'\n')
fss.close()

fsa=open(sys.argv[1]+'.sa', 'w')
count = 0
for line in open(sys.argv[1]+'.psa').readlines():
  count=count+1
  if count>13:
     cols = string.split(line)
     if (cols[2] in AA.keys()): fsa.write(cols[-3]+'\n')
fsa.close()

count = 0
fss=open(sys.argv[2]+'.ss', 'w')
for line in open(sys.argv[2]+'.ssm').readlines():
  count=count+1
  if count>13:
     cols = string.split(line)
     if (cols[2] in AA.keys()): fss.write(cols[-1]+'\n')
fss.close()

fsa=open(sys.argv[2]+'.sa', 'w')
count = 0
for line in open(sys.argv[2]+'.psa').readlines():
  count=count+1
  if count>13:
     cols = string.split(line)
     if (cols[2] in AA.keys()): fsa.write(cols[-3]+'\n')
fsa.close()
