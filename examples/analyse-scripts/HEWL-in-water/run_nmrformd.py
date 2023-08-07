import os
import numpy as np
import nmrformd as nmrmd
import MDAnalysis as mda

#import matplotlib.pyplot as plt

# Path to data
datapath = "../../raw-data/HEWL-in-water/"

alpha_m = [np.sqrt(16 * np.pi / 5), np.sqrt(8 * np.pi / 15), np.sqrt(32 * np.pi / 15)]

all_folders = ["T300K_ratio0.11/", "T300K_ratio0.61/"] #, "T300K_ratio2.73/"]

for folder in all_folders:

    # Import universe
    u = mda.Universe(datapath+folder+"conf.gro", datapath+folder+"prod.xtc")
    #u.transfer_to_memory(step=5, stop=50000)
    
    # Water
    water = u.select_atoms("name OW HW1 HW2")
    h_water = u.select_atoms("name HW1 HW2")

    # Protein
    all_name = ' '
    for name in np.unique(u.atoms.names):
        if (name != 'OW') & (name != 'HW2') & (name != 'HW1'):
            all_name += name + ' '
    protein = u.select_atoms('name '+all_name)
    all_name_H = ' '
    for name in np.unique(u.atoms.names):
        if (name != 'OW') & (name != 'HW2') & (name != 'HW1') & (name[0] == 'H'):
            all_name_H += name + ' '
    h_protein = u.select_atoms('name '+all_name_H)

    # All
    h_all = h_protein+h_water

    # Calculate NMR properties
    nmr_all_all = nmrmd.NMR(u, atom_group=h_all, neighbor_group=h_all, number_i=100) #np.min([40, h_all.atoms.n_atoms]))

    dictionary = {
        "nmr_all_f": nmr_all_all.f,
        "nmr_all_R1": nmr_all_all.R1,
        "nmr_all_R2": nmr_all_all.R2,
    }
    np.save(folder[:-1]+"_analysed-all-data.npy", dictionary)
