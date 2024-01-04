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


# In[4]:


while 1<2:
    # varying number
    mystep = 1
    for N in ["N25", "N39", "N62", "N99", "N158", "N251", "N398", "N631", "N1002", "N1589", "N2521"]:
        datapath = git_path + "/nmrformd-data/bulk-water/raw-data/"+N+"/"
        for n in np.arange(1, 10):
            if (os.path.exists(datapath+"prod"+str(n)+".xtc")) & (os.path.exists(datapath+"topology.data")):
                u = mda.Universe(datapath+"topology.data", datapath+"prod"+str(n)+".xtc")
                u.transfer_to_memory(step = mystep)
                hydrogen = u.select_atoms("type 2")
            elif (os.path.exists(datapath+"prod"+str(n)+".xtc")) & (os.path.exists(datapath+"prod1.tpr")):
                u = mda.Universe(datapath+"prod1.tpr", datapath+"prod"+str(n)+".xtc")
                u.transfer_to_memory(step = mystep)
                hydrogen = u.select_atoms("name HW1 HW2")
            dt = np.round(u.trajectory.dt,2)
            n_hydrogen = hydrogen.n_atoms
            intra = nmrformd.NMR(u, atom_group = hydrogen, neighbor_group = hydrogen, type_analysis = "intra_molecular", number_i=5)
            inter = nmrformd.NMR(u, atom_group = hydrogen, neighbor_group = hydrogen, type_analysis = "inter_molecular", number_i=1)
            n_intra = save_result(intra, name = N + "_intra_dt" + str(dt))
            n_inter = save_result(inter, name = N + "_inter_dt" + str(dt))
            if n == 1:
                print(n_hydrogen, mystep, n_intra)
    # varying step
    N = "N4000"
    for mystep in [1, 2, 4, 8, 16, 32]:
        datapath = git_path + "/nmrformd-data/bulk-water/raw-data/"+N+"/"
        xtcs = [[datapath+"prod1.xtc", datapath+"prod2.xtc"],
                [datapath+"prod3.xtc", datapath+"prod4.xtc"]]
        for xtc in xtcs:
            u = mda.Universe(datapath+"prod1.tpr", xtc)
            u.transfer_to_memory(step = mystep)
            dt = np.round(u.trajectory.dt,2)
            hydrogen = u.select_atoms("name HW1 HW2")
            n_hydrogen = hydrogen.n_atoms
            intra = nmrformd.NMR(u, atom_group = hydrogen, neighbor_group = hydrogen, type_analysis = "intra_molecular", number_i=5)
            inter = nmrformd.NMR(u, atom_group = hydrogen, neighbor_group = hydrogen, type_analysis = "inter_molecular", number_i=1)
            n_intra = save_result(intra, name = N + "_intra_dt" + str(dt))
            n_inter = save_result(inter, name = N + "_inter_dt" + str(dt))
            if mystep == 1:
                print(n_hydrogen, mystep, n_intra, dt)
    # high res
    N = "N4000_HR"
    for mystep in [1, 2, 4, 8, 16, 32]:
        datapath = git_path + "/nmrformd-data/bulk-water/raw-data/"+N+"/"
        xtc = [datapath+"prod1.xtc", datapath+"prod2.xtc",
               datapath+"prod3.xtc", datapath+"prod4.xtc",
               datapath+"prod5.xtc", datapath+"prod6.xtc",
               datapath+"prod7.xtc", datapath+"prod8.xtc",
               datapath+"prod9.xtc", datapath+"prod10.xtc"]
        u = mda.Universe(datapath+"prod1.tpr", xtc)
        u.transfer_to_memory(step = mystep)
        dt = np.round(u.trajectory.dt,2)
        hydrogen = u.select_atoms("name HW1 HW2")
        n_hydrogen = hydrogen.n_atoms
        intra = nmrformd.NMR(u, atom_group = hydrogen, neighbor_group = hydrogen, type_analysis = "intra_molecular", number_i=5)
        inter = nmrformd.NMR(u, atom_group = hydrogen, neighbor_group = hydrogen, type_analysis = "inter_molecular", number_i=1)
        n_intra = save_result(intra, name = N + "_intra_dt" + str(dt))
        n_inter = save_result(inter, name = N + "_inter_dt" + str(dt))
        if mystep == 1:
            print(n_hydrogen, mystep, n_intra, dt)
    # iso
    N = "N4000"
    datapath = git_path + "/nmrformd-data/bulk-water/raw-data/"+N+"/"
    xtcs = [[datapath+"prod1.xtc", datapath+"prod2.xtc"],
            [datapath+"prod3.xtc", datapath+"prod4.xtc"]]
    for xtc in xtcs:
        u = mda.Universe(datapath+"prod1.tpr", xtc)
        dt = np.round(u.trajectory.dt,2)
        hydrogen = u.select_atoms("name HW1 HW2")
        n_hydrogen = hydrogen.n_atoms
        intra = nmrformd.NMR(u, atom_group = hydrogen, neighbor_group = hydrogen,
                             type_analysis = "intra_molecular", number_i=5, isotropic = False)
        inter = nmrformd.NMR(u, atom_group = hydrogen, neighbor_group = hydrogen,
                             type_analysis = "inter_molecular", number_i=1, isotropic = False)
        n_intra = save_result(intra, name = N + "_intra_aniso_dt" + str(dt), aniso = True)
        n_inter = save_result(inter, name = N + "_inter_aniso_dt" + str(dt), aniso = True)
        print("N:", n_intra)

