import numpy as np
import os

def detect_groups(u):
    """Detects the hydrogens belonging to the water
    and to the protein. Returns the group."""
    # Water
    h_water = u.select_atoms("name HW1 HW2")
    # Protein
    all_name_H = ' '
    for name in np.unique(u.atoms.names):
        if (name != 'OW') & (name != 'HW2') & (name != 'HW1') & (name[0] == 'H'):
            all_name_H += name + ' '
    h_protein = u.select_atoms('name '+all_name_H)
    # All
    h_all = h_protein+h_water
    return h_all, h_water, h_protein

def save_nmr_data(filename, f_cut_off, nmr_h2o_intra, nmr_h2o_inter, nmr_pro_intra, nmr_h2o_pro, nmr_tot_tot):
    """Save the correlation functions in a dictionary"""
    N = 1
    f = nmr_h2o_intra.f[nmr_h2o_intra.f<f_cut_off]
    nmr_h2o_intra_R1 = nmr_h2o_intra.R1[nmr_h2o_intra.f<f_cut_off]
    nmr_h2o_intra_R2 = nmr_h2o_intra.R2[nmr_h2o_intra.f<f_cut_off]
    nmr_h2o_inter_R1 = nmr_h2o_inter.R1[nmr_h2o_intra.f<f_cut_off]
    nmr_h2o_inter_R2 = nmr_h2o_inter.R2[nmr_h2o_intra.f<f_cut_off]
    nmr_pro_intra_R1 = nmr_pro_intra.R1[nmr_h2o_intra.f<f_cut_off]
    nmr_pro_intra_R2 = nmr_pro_intra.R2[nmr_h2o_intra.f<f_cut_off]
    nmr_h2o_pro_R1 = nmr_h2o_pro.R1[nmr_h2o_intra.f<f_cut_off]
    nmr_h2o_pro_R2 = nmr_h2o_pro.R2[nmr_h2o_intra.f<f_cut_off]
    nmr_tot_tot_R1 = nmr_tot_tot.R1[nmr_h2o_intra.f<f_cut_off]
    nmr_tot_tot_R2 = nmr_tot_tot.R2[nmr_h2o_intra.f<f_cut_off]
    # if the data file exists already, merge all data
    if os.path.exists(filename):
        previous_dictionary = np.load(filename, allow_pickle=True)
        prev_nmr_h2o_intra_R1 = previous_dictionary.item()["nmr_h2o_intra_R1"]
        prev_nmr_h2o_intra_R2 = previous_dictionary.item()["nmr_h2o_intra_R2"]
        prev_nmr_h2o_inter_R1 = previous_dictionary.item()["nmr_h2o_inter_R1"]
        prev_nmr_h2o_inter_R2 = previous_dictionary.item()["nmr_h2o_inter_R2"]
        prev_nmr_pro_intra_R1 = previous_dictionary.item()["nmr_pro_intra_R1"]
        prev_nmr_pro_intra_R2 = previous_dictionary.item()["nmr_pro_intra_R2"]
        prev_nmr_h2o_pro_R1 = previous_dictionary.item()["nmr_h2o_pro_R1"]
        prev_nmr_h2o_pro_R2 = previous_dictionary.item()["nmr_h2o_pro_R2"]
        prev_nmr_tot_tot_R1 = previous_dictionary.item()["nmr_tot_tot_R1"]
        prev_nmr_tot_tot_R2 = previous_dictionary.item()["nmr_tot_tot_R2"]
        N_prev = np.real(previous_dictionary.item()["N"])
        nmr_h2o_intra_R1 = (nmr_h2o_intra_R1*N + prev_nmr_h2o_intra_R1*N_prev) / (N_prev + N)
        nmr_h2o_intra_R2 = (nmr_h2o_intra_R2*N + prev_nmr_h2o_intra_R2*N_prev) / (N_prev + N)
        nmr_h2o_inter_R1 = (nmr_h2o_inter_R1*N + prev_nmr_h2o_inter_R1*N_prev) / (N_prev + N)
        nmr_h2o_inter_R2 = (nmr_h2o_inter_R2*N + prev_nmr_h2o_inter_R2*N_prev) / (N_prev + N)
        nmr_pro_intra_R1 = (nmr_pro_intra_R1*N + prev_nmr_pro_intra_R1*N_prev) / (N_prev + N)
        nmr_pro_intra_R2 = (nmr_pro_intra_R2*N + prev_nmr_pro_intra_R2*N_prev) / (N_prev + N)
        nmr_h2o_pro_R1 = (nmr_h2o_pro_R1*N + prev_nmr_h2o_pro_R1*N_prev) / (N_prev + N)
        nmr_h2o_pro_R2 = (nmr_h2o_pro_R2*N + prev_nmr_h2o_pro_R2*N_prev) / (N_prev + N)
        nmr_tot_tot_R1 = (nmr_tot_tot_R1*N + prev_nmr_tot_tot_R1*N_prev) / (N_prev + N)
        nmr_tot_tot_R2 = (nmr_tot_tot_R2*N + prev_nmr_tot_tot_R2*N_prev) / (N_prev + N)
        N = N_prev + N
    # save data
    dictionary = {
        "f": f,
        "nmr_h2o_intra_R1": nmr_h2o_intra_R1,
        "nmr_h2o_intra_R2": nmr_h2o_intra_R2,
        "nmr_h2o_inter_R1": nmr_h2o_inter_R1,
        "nmr_h2o_inter_R2": nmr_h2o_inter_R2,
        "nmr_pro_intra_R1": nmr_pro_intra_R1,
        "nmr_pro_intra_R2": nmr_pro_intra_R2,
        "nmr_h2o_pro_R1": nmr_h2o_pro_R1,
        "nmr_h2o_pro_R2": nmr_h2o_pro_R2,
        "nmr_tot_tot_R1": nmr_tot_tot_R1,
        "nmr_tot_tot_R2": nmr_tot_tot_R2,
        "N": N,
    }
    np.save(filename, dictionary)