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
    u = mda.Universe("bulk_water_lammps/conf.data", "bulk_water_lammps/traj.xtc")
    group_i = "type 2"

    nmr_1 = NMR.NMR(u, [group_i, group_i],
                    type_analysis="full",
                    number_i = 40)
    assert nmr_1.T1 > 1
    assert nmr_1.T1 < 4

    nmr_2 = NMR.NMR(u, [group_i, group_i],
                    type_analysis="inter_molecular",
                    number_i = 40)
    assert nmr_2.T1 > 5.5
    assert nmr_2.T1 < 7.0

    nmr_3 = NMR.NMR(u, [group_i, group_i],
                    type_analysis="inter_molecular",
                    number_i = 40, f0 = 1e5)

    assert nmr_3.T1 > 36
    assert nmr_3.T1 < 38
