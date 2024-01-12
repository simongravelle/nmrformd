#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import MDAnalysis as mda
import nmrformd


# In[2]:


import sys, os, git
current_path = os.getcwd()
git_repo = git.Repo(current_path, search_parent_directories=True)
git_path = git_repo.git.rev_parse("--show-toplevel")


# In[3]:


def save_result(data, name = "intra_H2O"):
    """Save the correlation functions in dictionary"""
    if not os.path.exists("raw_data_long_time/"):
        os.makedirs("raw_data_long_time/")
    saving_file = "raw_data_long_time/" + name + ".npy"
    t = data.t
    f = data.f
    C = data.gij[0]
    R1 = data.R1
    R2 = data.R2
    N = data.group_j.atoms.n_atoms
    try:
        previous_dictionary = np.load(saving_file, allow_pickle=True)
        t_prev = np.real(previous_dictionary.item()["t"])
        assert len(t_prev) == len(t)
        C_prev = np.real(previous_dictionary.item()["C"])
        R1_prev = np.real(previous_dictionary.item()["R1"])
        R2_prev = np.real(previous_dictionary.item()["R2"])
        N_prev = np.real(previous_dictionary.item()["N"])
        C = (C*N + C_prev*N_prev) / (N_prev + N)
        R1 = (R1*N + R1_prev*N_prev) / (N_prev + N)
        R2 = (R2*N + R2_prev*N_prev) / (N_prev + N)
        N += N_prev
    except:
        pass
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


# In[4]:


folder = "/data/nmrformd-data/HEWL-in-water/raw-data/water-to-protein-0.73/"
tpr = folder+"prod.tpr"
all_xtc = [folder+"prod.xtc"]


# In[6]:


u = mda.Universe(tpr, all_xtc)
u.transfer_to_memory(step = 100)
duration = np.round(u.trajectory.dt * u.trajectory.n_frames,1) # ps
min_freq = 1/duration*1e6
max_freq = 1/(u.trajectory.dt)*1e6
print("minimum frequency", np.round(min_freq), "-- maximum frequency", np.round(max_freq))


# In[ ]:


# isolate the hydrogen from water and from the protein
H2O = u.select_atoms("name HW1 HW2")
HEWL = u.select_atoms("")
for name in np.unique(u.atoms.names):
    if (name[0] == "H") & (name != "HW1") & (name != "HW2"):
        HEWL += u.select_atoms("name "+name)
ALL = H2O + HEWL
while 1 < 2:
    # H2O-H2O
    for i in range(5):
        intra_H2O = nmrformd.NMR(u, atom_group=H2O, neighbor_group=H2O, type_analysis="intra_molecular", number_i=1)
        N = save_result(intra_H2O, name = "intra_H2O")
    inter_H2O = nmrformd.NMR(u, atom_group=H2O, neighbor_group=H2O, type_analysis="inter_molecular", number_i=1)
    full_H2O = nmrformd.NMR(u, atom_group=H2O, neighbor_group=H2O, number_i=1)
    print("Number of cycle water:", N)
    N = save_result(inter_H2O, name = "inter_H2O")
    N = save_result(full_H2O, name = "full_H2O")
    # HEWL-HEWL
    intra_HEWL = nmrformd.NMR(u, atom_group=HEWL, neighbor_group=HEWL, type_analysis="intra_molecular", number_i=1)
    inter_HEWL = nmrformd.NMR(u, atom_group=HEWL, neighbor_group=HEWL, type_analysis="inter_molecular", number_i=1)
    full_HEWL = nmrformd.NMR(u, atom_group=HEWL, neighbor_group=HEWL, number_i=1)
    N = save_result(intra_HEWL, name = "intra_HEWL")
    print("Number of cycle HEWL:", N)
    N = save_result(inter_HEWL, name = "inter_HEWL")
    N = save_result(full_HEWL, name = "full_HEWL")
    # H2O-HEWL
    inter_H2O_HEWL = nmrformd.NMR(u, atom_group=H2O, neighbor_group=HEWL, number_i=1)
    inter_HEWL_H2O = nmrformd.NMR(u, atom_group=HEWL, neighbor_group=H2O, number_i=1)
    # FULL
    FULL = nmrformd.NMR(u, atom_group=ALL, neighbor_group=ALL, number_i=1)
    # save_results
    N = save_result(inter_H2O_HEWL, name = "inter_H2O_HEWL")
    print("Number of cycle INTER:", N)
    N = save_result(inter_HEWL_H2O, name = "inter_HEWL_H2O")
    N = save_result(FULL, name = "FULL")

