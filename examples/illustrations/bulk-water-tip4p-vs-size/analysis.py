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

path = "/data/nmrformd-data/bulk-water-tip4p-vs-size/raw-data/"
run = True
R1_intra = 1e-5
R1_inter = 1e-5
while run:
    for N in [25, 39, 62, 99, 158, 251, 398, 631, 1002, 1589, 2521]:
        for n in range(10):
            # Import MDA universe
            datapath = path + "N"+str(N)+"/"
            if os.path.exists(datapath+"prod"+str(n)+".xtc"):
                # Read current status
                name_intra = "N"+str(N)+"_intra_n"+str(n)
                name_inter = "N"+str(N)+"_inter_n"+str(n) 
                n_intra = measure_N(name_intra)
                n_inter = measure_N(name_inter)
                u = mda.Universe(datapath+"topology.data", datapath+"prod"+str(n)+".xtc")
                hydrogen = u.select_atoms("type 2")
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
                    print("N", N, "no calculation")
                    #run = False
                # Print information
                R1_intra = measure_R1(name_intra)
                R1_inter = measure_R1(name_inter)
                print("N =", N, "n", n, "----", n_intra, n_hydrogen, 
                    "R1:", np.round(R1_intra,2), np.round(R1_inter,2), 
                    "T1:", np.round(1/(R1_intra+R1_inter),2))
    print(" ")
