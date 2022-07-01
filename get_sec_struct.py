import sys

from modeller import *
# from modeller.automodel import *

log.verbose()
env = environ()

#Use Modeller to get secondary structure and sovel accessibility values for the first and second structures

modl = model(env, file=(sys.argv[1]))
modl.write_data(file=(sys.argv[1]), output='SSM')
modl.write_data(file=(sys.argv[1]), output='PSA')

modl = model(env, file=(sys.argv[2]))
modl.write_data(file=(sys.argv[2]), output='SSM')
modl.write_data(file=(sys.argv[2]), output='PSA')


