#!/usr/bin/env python3
"""Main file for NMRforMD package."""
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2023 Authors and contributors
# Simon Gravelle
#
# Released under the GNU Public Licence, v3 or any higher version
# SPDX-License-Identifier: GPL-3.0-or-later

import random

import MDAnalysis.core.groups
import numpy as np
from scipy import constants as cst
from scipy.interpolate import interp1d
from scipy.special import sph_harm

from .utilities import autocorrelation_function, find_nearest, fourier_transform


class NMR:
    """Calculate NMR relaxation time from MDAnalysis universe.

    Parameters
    ----------
    u : MDAnalysis.Universe
        MDAnalysis universe containing all the information describing
        the molecular dynamics system.
    atom_group : MDAnalysis.AtomGroup
        Target atom groups for NMR calculation.
    neighbor_group : MDAnalysis.AtomGroup
        Neighbor atom groups. If not specified, atom_group is used.
    type_analysis : str, default ``full``
        Type of analysis, which can be ``full``, ``intra_molecular``,
        or ``inter_molecular``.
    number_i : int, default 0
        Number of atom of the target group to consider for the calculation.
        If ``number_i = 0``, all atoms are considered.
    isotropic : bool, default ``True``
        If isotropic is true, only the spherical harmonic of order 0 is considered, 
        which is usually valid for bulk systems. For non-isotropic systems,
        use ``False``.
    f0 : int, default ``None``
        Frequency at which ``T1`` and ``T2`` are calculated.
        If ``None``, ``f = 0`` is used.
    actual_dt : float, default ``None``
        Can be used to specify a different time interface between frames than the 
        one detected by MDAnalysis.
    hydrogen_per_atom : float, default 1.0
        Specify the number of hydrogen per atom, usefull for 
        coarse-grained simulations.
    pdb : bool, default True
        To turn off/on the periodic boundary condition treatment.            
    """

    def __init__(self,
                 u: MDAnalysis.Universe,
                 atom_group: MDAnalysis.AtomGroup,
                 neighbor_group: MDAnalysis.AtomGroup = None,
                 type_analysis: str = "full",
                 number_i: int = 0,
                 isotropic: bool = True,
                 f0: float = None,
                 actual_dt: float = None,
                 hydrogen_per_atom: float = 1.0,
                 spin: float = 1/2,
                 pbc: bool = True
                 ):
        
        """Initialise class NMR."""
        self.u = u
        self.target_i = atom_group
        if neighbor_group is None:
            self.neighbor_j = atom_group
        else:
            self.neighbor_j = neighbor_group
        self.type_analysis = type_analysis
        self.number_i = number_i
        self.f0 = f0
        self.isotropic = isotropic
        if self.isotropic:
            self.dim = 1
        else:
            self.dim = 3
        self.actual_dt = actual_dt
        self.hydrogen_per_atom = hydrogen_per_atom
        self.spin = spin
        self.pbc = pbc
        self.data = None
        self.rij = None
        self.r = None
        self.theta = None
        self.phi = None
        self.gij = None
        self.T1 = None
        self.T2 = None

        self.initialize()
        self.collect_data()
        self.finalize()

    def initialize(self):
        """Prepare the calculation"""
        self.verify_entry()
        self.define_constants()
        self.select_target_i()

    def verify_entry(self):
        """Verify that entries are correct, and that groups are not empty."""
        possible_analysis = ["inter_molecular", "intra_molecular", "full"]
        if self.type_analysis not in possible_analysis:
            raise ValueError("type_analysis can only be inter_molecular, intra_molecular, and full.")      
        if self.target_i.n_atoms == 0:
            raise ValueError("Empty atom group")
        if self.neighbor_j.n_atoms == 0:
            raise ValueError("Empty neighbor atom group")

    def define_constants(self):
        """Define prefactors.

        See this page for details : https://nmrformd.readthedocs.io
        Gamma is the gyromagnetic constant of 1H in Hz/T or C/kg
        K has the units of m^6/s^2
        alpha_m are normalizing coefficient for harmonic function
        """
        self.GAMMA = 2 * np.pi * 42.6e6 #todo offer this a an input parameter
        self.K = (3 / 2) * (cst.mu_0 / 4 / np.pi) ** 2 \
            * cst.hbar ** 2 * self.GAMMA ** 4 * self.spin * (1 + self.spin)
        self.alpha_m = [np.sqrt(16 * np.pi / 5),
                        np.sqrt(8 * np.pi / 15),
                        np.sqrt(32 * np.pi / 15)]

    def select_target_i(self):
        """Select the target atoms i
        
        If number_i=0, select all atoms in target_i
        If number_i > target_i.atoms.n_atoms, raise message
        Else, if 0<number_i<target_i.atoms.n_atoms, select atoms randomly
        """
        if self.number_i == 0:
            self.index_i = np.array(self.target_i.atoms.indices)
        elif self.number_i > self.target_i.atoms.n_atoms:
            print('Note : number_i is larger than the number of atoms in group target i\n'
                  '-> The number_i value will be ignored'
                  '-> All the atoms of the group i have been selected')
            self.index_i = np.array(self.target_i.atoms.indices)
        else:
            self.index_i = np.array(random.choices(self.target_i.atoms.indices, k=self.number_i))

    def collect_data(self):
        """Collect data by looping over atoms, time, and evaluate correlation"""
        # Loop on all the atom of group i
        for cpt_i, _ in enumerate(self.index_i):
            self.cpt_i = cpt_i
            #self.atom_index_i = atom_index_i # not used
            self.select_atoms_group_i()
            self.select_atoms_group_j()
            if cpt_i == 0:
                self.initialise_data()
            self.loop_over_trajectory()
            self.calculate_correlation_ij()

    def select_atoms_group_i(self):
        """Select atoms of the group i for the calculation."""
        self.group_i = self.u.select_atoms('index ' + str(self.index_i[self.cpt_i]))
        self.resids_i = self.group_i.resids[self.group_i.atoms.indices == self.index_i[self.cpt_i]]

    def select_atoms_group_j(self):
        """Select atoms of the group j for the calculation.

        For intra molecular analysis, group j are made of atoms of the
        same residue as group i.
        For inter molecular analysis, group j are made of atoms of
        different residues as group i.
        For full analysis, group j are made of atoms that are not in group i.
        """
        if self.type_analysis == "intra_molecular":
            same_residue : bool = self.neighbor_j.resids == self.resids_i
            different_atom : bool = self.neighbor_j.indices != self.index_i[self.cpt_i]
            index_j = self.neighbor_j.atoms.indices[same_residue & different_atom]
            str_j = ' '.join(str(e) for e in index_j)
        elif self.type_analysis == "inter_molecular":
            different_residue : bool = self.neighbor_j.resids != self.resids_i
            index_j = self.neighbor_j.atoms.indices[different_residue]
            str_j = ' '.join(str(e) for e in index_j)
        elif self.type_analysis == "full":
            different_atom : bool = self.neighbor_j.indices != self.index_i[self.cpt_i]
            index_j = self.neighbor_j.atoms.indices[different_atom]
            str_j = ' '.join(str(e) for e in index_j)
        if len(str_j) == 0:
            raise ValueError("Empty atom groups j \n"
                             "Wrong combination of type_analysis and group selection?")
        else:
            self.group_j = self.u.select_atoms('index ' + str_j)

    def initialise_data(self):
        """Initialise arrays.

        Create an array of zeros for the data and the correlation function.
        If anisotropic, the spherical harmonic may be complex, so dtype=complex64
        is used.
        Create an array for of values separated by timestep for the time. 
        """
        if self.isotropic:
            self.data = np.zeros((self.dim, self.u.trajectory.n_frames,
                                    self.group_j.atoms.n_atoms),
                                    dtype=np.float16)
        else:
            self.data = np.zeros((self.dim, self.u.trajectory.n_frames,
                                    self.group_j.atoms.n_atoms),
                                    dtype=np.complex64)
        self.gij = np.zeros((self.dim,  self.u.trajectory.n_frames),
                                dtype=np.float32)
        if self.actual_dt is None:
            self.timestep = np.round(self.u.trajectory.dt, 4)
        else:
            self.timestep = self.actual_dt
        self.t = np.arange(self.u.trajectory.n_frames) * self.timestep

    def loop_over_trajectory(self):
        """Loop of the MDA trajectory and extract rij. 
        
        Run over the MDA trajectory. If start, stop, or step are
        specified, only a sub-part of the trajectory is analyzed.
        """
        for cpt, ts in enumerate(self.u.trajectory):
            self.position_i = self.group_i.atoms.positions
            self.position_j = self.group_j.atoms.positions
            self.box = ts.dimensions
            # ensure that the box is orthonormal
            if not np.all(self.box[3:] == self.box[3:][0]) & np.all(self.box[3:][0] == 90.0):
                raise ValueError("NMRforMD does not accept non-orthogonal box"
                                 "You can use triclinic_to_orthorhombic from the package"
                                 "lipyphilic to convert the trajectory file.")
            self.vector_rij()
            self.cartesian_to_spherical()
            self.evaluate_function_F()
            self.data[:, cpt] = self.sph_val


    def vector_rij(self):
        """Calculate distance between position_i and position_j.
        
        By defaults, periodic boundary conditions are assumed. Pbc can be turned off using pbc = False.
        """
        if self.pbc:
            self.rij = (np.remainder(self.position_i - self.position_j
                         + self.box[:3]/2., self.box[:3]) - self.box[:3]/2.).T
        else:
            self.rij = (self.position_i - self.position_j).T

    def cartesian_to_spherical(self):
        """Convert cartesian coordinate to spherical."""
        self.r = np.sqrt(self.rij[0]**2 + self.rij[1]**2 + self.rij[2]**2)
        self.theta = np.arctan2(np.sqrt(self.rij[0]**2 + self.rij[1]**2), self.rij[2])
        self.phi = np.arctan2(self.rij[1], self.rij[0])

    def evaluate_function_F(self):
        """Evaluate the F functions.

        F = alpha Y / r ** 3
        Y : spherical harmonic
        r : spin-spin distance
        
        convention : theta = polar angle, phi = azimuthal angle
        note: scipy uses the opposite convention

        F has the units of Angstrom^(-6)
        """
        F_val = []
        for m in range(self.dim):
            F_val.append(self.alpha_m[m] * sph_harm(m, 2, self.phi, self.theta) / np.power(self.r, 3))
        if self.isotropic:
            F_val[0] = F_val[0].real
        self.sph_val = F_val

    def calculate_correlation_ij(self):
        """Calculate the correlation function."""
        for idx_j in range(self.group_j.atoms.n_atoms):
            for m in range(self.dim):
                self.gij[m] += autocorrelation_function(self.data[m, :, idx_j])

        self.gij = np.real(self.gij)

    def finalize(self):
        # calculate spectrums
        self.normalize_Gij()
        self.calculate_fourier_transform()
        self.calculate_spectrum()
        self.calculate_relaxationtime()

    def normalize_Gij(self):
        """Divide Gij by the number of spin pairs.
        
        Optional, for coarse grained model, apply a coefficient "hydrogen_per_atom" != 1
        """
        # normalise gij by the number of iteration (or number of pair spin)
        self.gij /= self.cpt_i+1
        if self.hydrogen_per_atom != 1:
            self.gij *= np.float32(self.hydrogen_per_atom)

    def calculate_fourier_transform(self):
        """Calculate spectral density J.
        
        Calculate the spectral density J from the 
        Fourier transform of the correlation function.
        """
        # for coarse grained models, possibly more than 1 hydrogen per atom
        self.J = []
        for m in range(self.dim):
            fij = fourier_transform(np.vstack([self.t, self.gij[m]]).T)
            self.J.append(np.real(fij.T[1]))
        self.J = np.array(self.J)
        self.f = np.real(fij.T[0])

    def calculate_spectrum(self):
        """Calculate spectrums R1 and R2 from J."""
        inter1d_0 = interp1d(self.f, self.J[0], fill_value="extrapolate")
        if self.isotropic:
            self.R1 = self.K/cst.angstrom ** 6 * (inter1d_0(self.f)
                                                  + 4 * inter1d_0(2 * self.f) )/6
            self.R2 = self.K/cst.angstrom ** 6 * (3/2*inter1d_0(self.f[0])
                                                  + 5/2*inter1d_0(self.f)
                                                  + inter1d_0(2 * self.f) )/6
        elif self.isotropic is False:
            inter1d_1 = interp1d(self.f, self.J[1], fill_value="extrapolate")
            inter1d_2 = interp1d(self.f, self.J[2], fill_value="extrapolate")
            self.R1 = self.K/cst.angstrom ** 6 * (inter1d_1(self.f)
                                                  + inter1d_2(2 * self.f))
            self.R2 = self.K/cst.angstrom ** 6 * ((1/4)*(inter1d_0(self.f[0])
                                                         + 10*inter1d_1(self.f)
                                                         + inter1d_2(2 * self.f)))
        
    def calculate_relaxationtime(self):
        """Calculate the relaxation time at a given frequency f0 (default is 0)"""
        if self.f0 is None:
            self.T1 = 1/self.R1[0]
            self.T2 = 1/self.R2[0]
        else:
            idx = find_nearest(self.f, self.f0)
            self.T1 = 1 / self.R1[idx]
            self.T2 = 1 / self.R2[idx]
