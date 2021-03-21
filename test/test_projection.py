import os
import json
import sys
sys.path.append(os.path.dirname(os.getcwd()))
from src.pyAOM_utils import *

with open('single_molecules.json') as fp:
    ref_data=json.load(fp)

# STO parameters for GTO-STO projection step
STO_proj_dict={
    'H':{'STOs':1  ,'1s':1.0000},
    'C':{'STOs':1+3,'2s':1.6083,'2p':1.442657},
    'N':{'STOs':1+3,'2s':1.9237,'2p':1.646703},
    'O':{'STOs':1+3,'2s':2.2458,'2p':1.858823},
    'F':{'STOs':1+3,'2s':2.5638,'2p':2.136394},
    'S':{'STOs':1+3,'3s':2.1223,'3p':1.651749},
}

test_data=projection_reg_test(ref_data,STO_proj_dict)