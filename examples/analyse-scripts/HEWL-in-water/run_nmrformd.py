#!/usr/bin/env python
# coding: utf-8

# Import libraries
import os
import numpy as np
import nmrformd as nmrmd
import MDAnalysis as mda

# Path to data
datapath = "../../raw-data/HEWL-in-water/"

alpha_m = [np.sqrt(16 * np.pi / 5), np.sqrt(8 * np.pi / 15), np.sqrt(32 * np.pi / 15)]

all_folders = ["T300K_ratio0.61/", "T300K_ratio2.73/"]

for folder in all_folders:

    # Import universe
    u = mda.Universe(datapath+folder+"conf.gro", datapath+folder+"prod.xtc")
    
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

    # Calculate NMR properties
    nmr_water = nmrmd.NMR(u, h_water, number_i=np.min([100, h_water.atoms.n_atoms]))
    print(folder, "water calculated")
    nmr_protein = nmrmd.NMR(u, h_protein, number_i=100)
    print(folder, "protein calculated")
    print()

    dictionary = {
        "nmr_water_f": nmr_water.f,
        "nmr_water_R1": nmr_water.R1,
        "nmr_water_R2": nmr_water.R2,
        "nmr_protein_f": nmr_protein.f,
        "nmr_protein_R1": nmr_protein.R1,
        "nmr_protein_R2": nmr_protein.R2,
    }
    np.save(datapath+folder+"analysed-data.npy", dictionary)

