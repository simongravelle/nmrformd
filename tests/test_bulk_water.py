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
    u = mda.Universe("../examples/bulk-water/lammps-inputs/topology.data",
                        "../examples/bulk-water/lammps-inputs/traj.xtc")

    group_i = u.select_atoms("type 2")
    nmr_1 = nmrmd.NMR(u, group_i, type_analysis="inter_molecular",
                      isotropic=True, number_i=10)
    assert np.isclose(nmr_1.T1 , nmr_1.T2)
    # expected values is near 6-7 seconds
    assert (nmr_1.T1 > 3) & (nmr_1.T1 < 10)