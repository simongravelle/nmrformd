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

import nmrformd as nmrmd


def test_nmr():
    """Test NMR module using bulk water trajectory."""
    u = mda.Universe("bulk_water_lammps/conf.data",
                     "bulk_water_lammps/traj.xtc")
    ids = ""
    for i in np.arange(1, 300, 24):
        ids = ids + " " + str(i)
    group_i = "type 2 and index"+ids
    group_j = "type 2"
    nmr_1 = nmrmd.NMR(u, [group_i, group_j],
                    type_analysis="inter_molecular")
    assert np.isclose(nmr_1.tau, 4.66710)
    assert np.isclose(nmr_1.delta_omega, 37.80930)
    assert np.isclose(nmr_1.T1, 6.833889)
    assert np.isclose(nmr_1.T2, 6.833889)

    nmr_2 = nmrmd.NMR(u, [group_i, group_j],
                    type_analysis="inter_molecular",
                    order="m012")
    assert np.isclose(nmr_2.tau[0], 4.667102)
    assert np.isclose(nmr_2.delta_omega[0], 37.809309)
    assert np.isclose(nmr_2.T1, 7.1193362)

    nmr_3 = nmrmd.NMR(u, [group_i, group_j], f0=1e5)
    assert np.isclose(nmr_3.T1, 11.42645)
    assert np.isclose(nmr_3.T2, 3.677240)
