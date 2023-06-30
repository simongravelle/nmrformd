#!/usr/bin/env python

import MDAnalysis as mda
import numpy as np
import nmrformd

def save_file(result, file_name, group1, group2):
    dictionary = {
    "t": result.t,
    "correlation_function": result.gij,
    "T1": result.T1,
    "group1": group1,             
    "group2": group2
    }
    np.save(file_name, dictionary) 

# import the trajectory
tpr = "run.tpr"
xtc = "run.xtc"
u = mda.Universe(tpr, xtc)

# create groups (only water molecules here)
group1 = "type HW"
group2 = "type HW HCA2 HCP1"
# uncomment to select PEG instead
# group1 = "type HCA2 HCP1"
# group2 = "type HW HCA2 HCP1"
# uncomment to select ALL instead
# group1 = "type HW HCA2 HCP1"
# group2 = "type HW HCA2 HCP1"
print("group 1 is " + group1)
print("group 2 is " + group2)

# run NMRforMD
# intermolecular analysis
result = nmrformd.NMR(u, [group1, group2], type_analysis="inter_molecular")
save_file(result, "H2O_intermolecular.npy", group1, group2)

# intramolecular analysis
result = nmrformd.NMR(u, [group1, group2], type_analysis="intra_molecular")
save_file(result, "H2O_intramolecular.npy", group1, group2)

