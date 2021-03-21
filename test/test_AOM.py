import os
import json
import sys
sys.path.append(os.path.dirname(os.getcwd()))
from src.pyAOM_utils import *
import time

with open('dimers.json') as fp:
    ref_data=json.load(fp)
    
# define STO parameters for AOM calculations
AOM_dict={
    'H':{'STOs':1  ,'1s':1.0000},
    'C':{'STOs':1+3,'2s':1.6083,'2p':1.385600},
    'N':{'STOs':1+3,'2s':1.9237,'2p':1.617102},
    'O':{'STOs':1+3,'2s':2.2458,'2p':1.505135},
    'F':{'STOs':1+3,'2s':2.5638,'2p':1.665190},
    'S':{'STOs':1+3,'3s':2.1223,'3p':1.641119},
}

# HAB79 regression test
test_data=overlap_reg_test(ref_data,AOM_dict)