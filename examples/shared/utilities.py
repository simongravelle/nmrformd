import numpy as np
import os

def save_result(data, name, aniso = False):
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

def measure_N(name):
    """measure_N"""
    try:
        saving_file = "raw_data/" + name + ".npy"
        previous_dictionary = np.load(saving_file, allow_pickle=True)
        N = np.real(previous_dictionary.item()["N"])
    except:
        N = 0
    return N

def measure_R1(name):
    """measure_N"""
    try:
        saving_file = "raw_data/" + name + ".npy"
        previous_dictionary = np.load(saving_file, allow_pickle=True)
        R1 = np.real(previous_dictionary.item()["R1"])
        return R1[0]
    except:
        return 1e-5