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

#from datafiles import (SPCE_GRO, SPCE_ITP)

#from pkg_resources import resource_filename

import os
cwd = os.getcwd()
print(cwd)

mda.Universe("water/spce.itp", "water/spce.gro", topology_format='itp')

def create_universe(n_molecules, angle_deg):
    """Create universe with regularly-spaced water molecules."""
    fluid = []
    for _n in range(n_molecules):
        fluid.append(mda.Universe("water/spce.itp",
                                  "water/spce.gro",
                                  topology_format='itp'))
    dimensions = fluid[0].dimensions

    translations = [(0, 0, 0),
                    (0, 0, 10)]

    for molecule, translation in zip(fluid, translations):
        molecule.atoms.translate(translation)
    u = mda.Merge(*[molecule.atoms for molecule in fluid])

    dimensions[2] *= n_molecules
    u.dimensions = dimensions
    u.residues.molnums = list(range(1, n_molecules + 1))
    return u

class TestNMR(object):
    """Tests for the NMR class using only two water molecules.

    The distance between the two molecules is 1 nm, and
    they do not move.
    """

    def test_2_water_0(self):
        u =  create_universe(2, 0)
        group_H2O_1 = ["name HW1"]

        #print(group_H2O_1.atoms.names)
        #print(group_H2O_1.atoms.positions)

        nmr_result = NMR.NMR(u, group_H2O_1, "full", 0, "m0", actual_dt = 2)
        