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
    """Test NMR module using bulk water trajectory."""
    u = mda.Universe("bulk_water_lammps/conf.data",
                     "bulk_water_lammps/traj.xtc")
    ids = ""
    for i in np.arange(1, 300, 6):
        ids = ids + " " + str(i)
    group_i = "type 2 and index"+ids
    group_j = "type 2"
    nmr_1 = NMR.NMR(u, [group_i, group_j],
                    type_analysis="inter_molecular")
    assert np.isclose(nmr_1.tau, 4.796872)
    assert np.isclose(nmr_1.delta_omega, 38.22099)
    assert np.isclose(nmr_1.T1, 6.50655)

    nmr_2 = NMR.NMR(u, [group_i, group_j],
                    type_analysis="inter_molecular",
                    order="m012")
    assert np.isclose(nmr_2.tau[0], 4.796872)
    assert np.isclose(nmr_2.delta_omega[0], 38.22099)
    assert np.isclose(nmr_2.T1, 6.99148)
