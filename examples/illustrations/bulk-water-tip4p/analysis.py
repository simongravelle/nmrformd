#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import MDAnalysis as mda
import nmrformd


# In[ ]:


import sys, os, git
current_path = os.getcwd()
git_repo = git.Repo(current_path, search_parent_directories=True)
git_path = git_repo.git.rev_parse("--show-toplevel")


# In[ ]:


def save_result(data, name = "intra_H2O", aniso = False):
    """Save the correlation functions in dictionary"""
    if not os.path.exists("raw_data/"):
        os.makedirs("raw_data/")
    saving_file = "raw_data/" + name + ".npy"
    t = data.t
    f = data.f
    if aniso:
        C1 = data.gij[0]
        C2 = data.gij[1]
        C3 = data.gij[2]
    else:
        C = data.gij[0]
    R1 = data.R1
    R2 = data.R2
    N = data.group_j.atoms.n_atoms
    try:
        previous_dictionary = np.load(saving_file, allow_pickle=True)
        t_prev = np.real(previous_dictionary.item()["t"])
        assert len(t_prev) == len(t)
        if aniso:
            C1_prev = np.real(previous_dictionary.item()["C1"])
            C2_prev = np.real(previous_dictionary.item()["C2"])
            C3_prev = np.real(previous_dictionary.item()["C3"])
        else:
            C_prev = np.real(previous_dictionary.item()["C"])
        R1_prev = np.real(previous_dictionary.item()["R1"])
        R2_prev = np.real(previous_dictionary.item()["R2"])
        N_prev = np.real(previous_dictionary.item()["N"])
        if aniso:
            C1 = (C1*N + C1_prev*N_prev) / (N_prev + N)
            C2 = (C2*N + C2_prev*N_prev) / (N_prev + N)
            C3 = (C3*N + C3_prev*N_prev) / (N_prev + N)
        else:
            C = (C*N + C_prev*N_prev) / (N_prev + N)
        R1 = (R1*N + R1_prev*N_prev) / (N_prev + N)
        R2 = (R2*N + R2_prev*N_prev) / (N_prev + N)
        N += N_prev
    except:
        pass
    if aniso:
        dictionary = {
        "t": t,
        "f": f,
        "C1": C1,
        "C2": C2,
        "C3": C3,
        "N": N,
        "R1": R1,
        "R2": R2,
        }
    else:
        dictionary = {
        "t": t,
        "f": f,
        "C": C,
        "N": N,
        "R1": R1,
        "R2": R2,
        }
    np.save(saving_file, dictionary)
    return N


# In[ ]:


# In[8]:


while 1<2:   
    # vary T
    for T in [280, 290, 310, 320]:
        datapath = git_path + "/data/bulk-water-tip4p/bulk-water/raw-data/N4000_" + str(T) + "K/"
        xtc = [datapath+"prod1.xtc", datapath+"prod2.xtc"]
        u = mda.Universe(datapath+"prod1.tpr", xtc)
        hydrogen = u.select_atoms("name HW1 HW2")
        n_hydrogen = hydrogen.n_atoms
        intra = nmrformd.NMR(u, atom_group = hydrogen, neighbor_group = hydrogen, type_analysis = "intra_molecular", number_i=5)
        inter = nmrformd.NMR(u, atom_group = hydrogen, neighbor_group = hydrogen, type_analysis = "inter_molecular", number_i=1)
        n_intra = save_result(intra, name = "N4000_intra_T" + str(T) + "K")
        n_inter = save_result(inter, name = "N4000_inter_T" + str(T) + "K")
        print(n_hydrogen, T, n_intra)
    # vary co
    for co in [6, 8, 10, 12]:
        datapath = git_path + "/data/bulk-water-tip4p/bulk-water/raw-data/N4000_co" + str(co) + "/"
        xtc = [datapath+"prod1.xtc", datapath+"prod2.xtc"]
        u = mda.Universe(datapath+"prod1.tpr", xtc)
        hydrogen = u.select_atoms("name HW1 HW2")
        n_hydrogen = hydrogen.n_atoms
        intra = nmrformd.NMR(u, atom_group = hydrogen, neighbor_group = hydrogen, type_analysis = "intra_molecular", number_i=5)
        inter = nmrformd.NMR(u, atom_group = hydrogen, neighbor_group = hydrogen, type_analysis = "inter_molecular", number_i=1)
        n_intra = save_result(intra, name = "N4000_intra_co" + str(co))
        n_inter = save_result(inter, name = "N4000_inter_co" + str(co))
        print(n_hydrogen, co, n_intra)

