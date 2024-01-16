#!/usr/bin/env python
# coding: utf-8

import numpy as np
import MDAnalysis as mda
import nmrformd

import sys, os, git
current_path = os.getcwd()
git_repo = git.Repo(current_path, search_parent_directories=True)
git_path = git_repo.git.rev_parse("--show-toplevel")

sys.path.append(git_path + '/examples/shared/')
from utilities import save_result, measure_N, measure_R1

path = "/data/nmrformd-data/bulk-water-tip4p-vs-temperature/raw_data/"
run = False
R1_intra = 1e-5
R1_inter = 1e-5
while run:
    for T in np.arange(280, 330, 10):
        # Read current status
        name_intra = "N4000_intra_T" + str(T) + "K"
        name_inter = "N4000_inter_T" + str(T) + "K"
        n_intra = measure_N(name_intra)
        n_inter = measure_N(name_inter)
        # Import MDA universe
        datapath = path + "N4000_" + str(T) + "K/"
        u = mda.Universe(datapath+"prod.tpr", datapath+"prod.xtc")
        hydrogen = u.select_atoms("name HW1 HW2")
        n_hydrogen = hydrogen.n_atoms
        if n_intra < n_hydrogen:
            # Calculated NMR properties
            intra = nmrformd.NMR(u, atom_group = hydrogen,
                                neighbor_group = hydrogen,
                                type_analysis = "intra_molecular",
                                number_i = 1)
            save_result(intra, name = name_intra)
            if np.random.random() < 0.2:
                inter = nmrformd.NMR(u, atom_group = hydrogen,
                                    neighbor_group = hydrogen,
                                    type_analysis = "inter_molecular",
                                    number_i = 1)
                save_result(inter, name = "N4000_inter_T" + str(T) + "K")
        else:
            run = False
        # Print information
        R1_intra = measure_R1(name_intra)
        R1_inter = measure_R1(name_inter)
        print("T =", T, "K --", n_intra, n_hydrogen, 
            "R1:", np.round(R1_intra,2), np.round(R1_inter,2), 
            "T1:", np.round(1/(R1_intra+R1_inter),2))
    print(" ")

# long trajectory
path = "/data/nmrformd-data/bulk-water-tip4p-vs-temperature/raw_data/"
run = True
R1_intra = 1e-5
R1_inter = 1e-5
while run:
    T = 300
    # Read current status
    name_intra = "N4000_intra_T" + str(T) + "K_long"
    name_inter = "N4000_inter_T" + str(T) + "K_long"
    n_intra = measure_N(name_intra)
    n_inter = measure_N(name_inter)
    # Import MDA universe
    datapath = path + "N4000_" + str(T) + "K_long/"
    u = mda.Universe(datapath+"prod.tpr", datapath+"prod.xtc")
    hydrogen = u.select_atoms("name HW1 HW2")
    n_hydrogen = hydrogen.n_atoms
    if n_intra < n_hydrogen:
        # Calculated NMR properties
        intra = nmrformd.NMR(u, atom_group = hydrogen,
                            neighbor_group = hydrogen,
                            type_analysis = "intra_molecular",
                            number_i = 1)
        save_result(intra, name = name_intra)
        if np.random.random() < 0.2:
            inter = nmrformd.NMR(u, atom_group = hydrogen,
                                neighbor_group = hydrogen,
                                type_analysis = "inter_molecular",
                                number_i = 1)
            save_result(inter, name = name_inter)
    else:
        run = False
    # Print information
    R1_intra = measure_R1(name_intra)
    R1_inter = measure_R1(name_inter)
    print("T =", T, "K --", n_intra, n_hydrogen, 
        "R1:", np.round(R1_intra,2), np.round(R1_inter,2), 
        "T1:", np.round(1/(R1_intra+R1_inter),2))
print(" ")
