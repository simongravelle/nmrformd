#!/usr/bin/env python3
"""Test file for NMRForMD package."""
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2022 Authors and contributors
# Simon Gravelle
#
# Released under the GNU Public Licence, v3 or any higher version
# SPDX-License-Identifier: GPL-3.0-or-later

import MDAnalysis as mda
import numpy as np

import nmrformd as NMR


def test_nmr():
    """Test NMR module."""
    u = mda.Universe("bulk_h2o/topology.tpr", "bulk_h2o/trajectory.xtc")

    group_i = "type HW  and index 0:10"
    group_j = "type HW and index 0:10"
    type_analysis = "intra_molecular"
    number_i = 0
    order = "m0"

    atom_group = [group_i, group_j]

    nmr_result = NMR.NMR(u, atom_group, type_analysis,
                         number_i, order)
    T1 = nmr_result.T1

    print("\n")
    print("\n")
    print("T1 = " + str(np.round(T1, 3)) + " s\n")
    print("\n")
    print("\n")


    #nmr_result = NMR.NMR(u, atom_group, "inter_molecular",
    #                     number_i, order)
    #T1 = nmr_result.T1
    #print("T1 = " + str(np.round(T1, 3)) + " s")


    #nmr_result = NMR.NMR(u, atom_group, "full",
    #                     number_i, order)
    #T1 = nmr_result.T1
    #print("T1 = " + str(np.round(T1, 3)) + " s")

    #assert np.isclose(T1, 1.029287074574358, rtol=1e-3, atol=0)



    #assert np.isclose(T1, 1.029287074574358, rtol=1e-3, atol=0)
    #assert np.isclose(T2, 1.029287074574358, rtol=1e-3, atol=0)
    #assert np.isclose(Dw, 56.996441578600304, rtol=1e-3, atol=0)
    #assert np.isclose(tau, 13.635804446316115, rtol=1e-3, atol=0)

    #group_i = "type HW and index 0:20"
    #group_j = "type HW and index 0:20"
    #analysis_type = "inter_molecular"
    #number_i = 0
    #order = "m0"
    #nmr_result = NMR.NMR(u, group_i,
    #                     group_j, analysis_type,
    #                     number_i, order)
    #T1 = nmr_result.T1
    #T2 = nmr_result.T2
    #Dw = nmr_result.delta_omega
    #tau = nmr_result.tau

    #assert np.isclose(T1, 1.0305, rtol=1e-4, atol=0)
    #assert np.isclose(T2, 1.0305, rtol=1e-4, atol=0)
    #assert np.isclose(Dw, 56.7908, rtol=1e-4, atol=0)
    #assert np.isclose(tau, 13.7180, rtol=1e-4, atol=0)

    #group_i = "type HW and index 0:40"
    #group_j = "type HW and index 0:40"
    #analysis_type = "intra_molecular"
    #number_i = 0
    #order = "m0"
    #nmr_result = NMR.NMR(u, group_i,
    #                     group_j, analysis_type,
    #                     number_i, order)
    #T1 = nmr_result.T1
    #T2 = nmr_result.T2
    #Dw = nmr_result.delta_omega
    #tau = nmr_result.tau

    #assert np.isclose(T1, 400.4619, rtol=1e-4, atol=0)
    #assert np.isclose(T2, 400.4619, rtol=1e-4, atol=0)
    #assert np.isclose(Dw, 5.99748, rtol=1e-4, atol=0)
    #assert np.isclose(tau, 3.16528, rtol=1e-4, atol=0)

    #group_i = "type HW and index 0:40"
    #group_j = "type HW and index 0:40"
    #analysis_type = "full"
    #number_i = 0
    #order = "m012"
    #f0 = 10000
    #nmr_result = NMR.NMR(u, group_i,
    #                     group_j, analysis_type,
    #                     number_i, order, f0)
    #T1 = nmr_result.T1
    #T2 = nmr_result.T2
    #Dw = nmr_result.delta_omega
    #tau = nmr_result.tau

    #assert np.isclose(T1, 4.931728, rtol=1e-4, atol=0)
    #assert np.isclose(T2, 3.853544, rtol=1e-4, atol=0)
    #assert np.isclose(Dw[0], 56.152647, rtol=1e-4, atol=0)
    #assert np.isclose(tau[0], 6.891118, rtol=1e-4, atol=0)

    #assert 1 == 3

test_nmr()