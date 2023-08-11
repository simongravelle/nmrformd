#!/usr/bin/env python
# coding: utf-8

# In[6]:


import os
import numpy as np
import nmrformd as nmrmd
import MDAnalysis as mda

import matplotlib.pyplot as plt

from utilities import detect_groups, save_nmr_data

# Path to data
datapath = "../../raw-data/HEWL-in-water/"

f_cut_off = 5e4

for rho_water in [0.08, 0.17, 0.25, 0.39, 0.52, 0.73]:
    folder = "T300K_ratio"+str(rho_water)+"/"
    filename = "T300K_ratio"+str(rho_water)+".npy"

    # Import universe
    u = mda.Universe(datapath+folder+"prod.tpr", datapath+folder+"prod.xtc")
    h_all, h_water, h_protein = detect_groups(u)

    # Calculate NMR properties for only one atom at a time
    nmr_h2o_intra = nmrmd.NMR(u, atom_group=h_water, neighbor_group=h_water, number_i=1, type_analysis="intra_molecular")
    nmr_h2o_inter = nmrmd.NMR(u, atom_group=h_water, neighbor_group=h_water, number_i=1, type_analysis="inter_molecular")
    nmr_pro_intra = nmrmd.NMR(u, atom_group=h_protein, neighbor_group=h_protein, number_i=1)
    nmr_h2o_pro = nmrmd.NMR(u, atom_group=h_water, neighbor_group=h_protein, number_i=1)
    nmr_tot_tot = nmrmd.NMR(u, atom_group=h_all, neighbor_group=h_all, number_i=1)

    save_nmr_data(filename, f_cut_off, nmr_h2o_intra, nmr_h2o_inter, nmr_pro_intra, nmr_h2o_pro, nmr_tot_tot)

