#!/usr/bin/env python3
"""Tests for the NMR class using only two water molecules."""
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2022 Authors and contributors
# Simon Gravelle
#
# Released under the GNU Public Licence, v3 or any higher version
# SPDX-License-Identifier: GPL-3.0-or-later

import MDAnalysis as mda
from scipy import constants as cst
import numpy as np
import nmrformd as nmrmd

def import_universe():
    data = "two_water/twomolecules.data"
    trj = "two_water/twomolecules.xtc"
    u = mda.Universe(data, trj)
    group_H2O = u.select_atoms("type 2")
    return u, group_H2O

def test_shape_array():
    """Assert that the data matrix dimension is consistent."""
    u, group_H2O = import_universe()
    n_frames = u.trajectory.n_frames
    nmr_1 = nmrmd.NMR(u, group_H2O, isotropic="True")
    assert nmr_1.data.shape == (1, n_frames, 1)
    assert nmr_1.gij.shape == (1, n_frames)
    assert nmr_1.t.shape == (n_frames,)
    n_fft = len(np.arange(np.int32(n_frames / 2)+1))
    assert nmr_1.J[0].shape == (n_fft,)
    assert nmr_1.f.shape == (n_fft,)
    assert nmr_1.R1.shape == (n_fft,)
    assert nmr_1.R2.shape == (n_fft,)
    nmr_2 = nmrmd.NMR(u, group_H2O, isotropic=False)
    assert nmr_2.data.shape == (3, n_frames, 1)
    assert nmr_2.gij.shape == (3, n_frames)
    assert nmr_2.t.shape == (n_frames,)
    assert nmr_2.J.shape == (3, n_fft)
    assert nmr_2.f.shape == (n_fft,)
    assert nmr_2.R1.shape == (n_fft,)
    assert nmr_2.R2.shape == (n_fft,)
    assert nmr_2.gij.T[0].shape == (3,)

def test_distance():
    """Assert that the calculated distance is correct."""
    u, group_H2O = import_universe()
    nmr_1 = nmrmd.NMR(u, group_H2O)
    assert np.isclose(nmr_1.rij[0][0], -8.0)
    assert np.isclose(nmr_1.rij[1][0], 0.0)
    assert np.isclose(nmr_1.rij[2][0], 0.0)
    nmr_2 = nmrmd.NMR(u, group_H2O, pbc = False)
    assert np.isclose(nmr_2.rij[0][0], 10.0)
    assert np.isclose(nmr_2.rij[1][0], 0.0)
    assert np.isclose(nmr_2.rij[2][0], 0.0)

def test_spherical():
    u, group_H2O = import_universe()
    nmr_1 = nmrmd.NMR(u, group_H2O)
    assert np.isclose(nmr_1.r[0], 8.0)
    assert np.isclose(nmr_1.theta[0], np.pi/2)
    assert np.isclose(nmr_1.phi[0], np.pi)
    nmr_2 = nmrmd.NMR(u, group_H2O, pbc = False)
    assert np.isclose(nmr_2.r[0], 10.0)
    assert np.isclose(nmr_2.theta[0], np.pi/2)
    assert np.isclose(nmr_2.phi[0], 0)

def test_correlation():
    u, group_H2O = import_universe()
    nmr_1 = nmrmd.NMR(u, group_H2O)
    #_pref = nmr_1._K / cst.angstrom ** 6
    #assert np.isclose(np.unique(nmr_1.gij/_pref), 1/8**6)
    #nmr_2 = nmrmd.NMR(u, group_H2O_1, order="m012")
    #assert np.isclose(np.unique(nmr_2.gij[0] / _pref), 1.0 / 8 ** 6)
    #assert np.isclose(np.unique(nmr_2.gij[1] / _pref), 0)
    #assert np.isclose(np.unique(nmr_2.gij[2] / _pref), 1.0 / 8 ** 6)
    #nmr_3 = nmrmd.NMR(u, group_H2O_1, pbc=False)
    #assert np.isclose(np.unique(nmr_3.gij / _pref), 1.0 / 10 ** 6)
    #nmr_4 = nmrmd.NMR(u, group_H2O_1, hydrogen_per_atom=2)
    #assert np.isclose(np.unique(nmr_4.gij / _pref), 2.0 / 8 ** 6)

def test_fourier_transform():
    u, group_H2O = import_universe()
    nmr_1 = nmrmd.NMR(u, group_H2O)
    #assert np.isclose(1e8 * np.unique(nmr_1.J[0]), 4.89951158)
    #assert np.isclose(1e8 * np.unique(nmr_1.J[1]), 0.0)

def test_group():
    """Assert that the number of atoms in both groups in correct."""
    u, group_H2O = import_universe()
    #nmr_1 = nmrmd.NMR(u, group_H2O, type_analysis="intra_molecular")
    #assert nmr_1.group_i.atoms.n_atoms == 1
    #assert nmr_1.group_j.atoms.n_atoms == 1
    #nmr_2 = nmrmd.NMR(u, group_H2O_1, type_analysis="inter_molecular")
    #assert nmr_2.group_i.atoms.n_atoms == 1
    #assert nmr_2.group_j.atoms.n_atoms == 2
    #nmr_3 = nmrmd.NMR(u, group_H2O_1, type_analysis="full")
    #assert nmr_3.group_i.atoms.n_atoms == 1
    #assert nmr_3.group_j.atoms.n_atoms == 3