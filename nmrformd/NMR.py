#!/usr/bin/env python3
"""Main file for NMRforMD package."""
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2022 Authors and contributors
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
from MDAnalysis.analysis.distances import distance_array

from .utilities import correlation_function, find_nearest, fourier_transform


class NMR:
    """Calculate NMR relaxation time from MDAnalysis universe.

    Parameters
    ----------
    u : MDAnalysis.Universe
        MDAnalysis universe containing all the information describing
        the molecular dynamics system.
    atom_group : [MDAnalysis.AtomGroup, MDAnalysis.AtomGroup]
        Target and Neighbor groups, respectively.
    type_analysis : str, default ``full``
        Type of analysis, which can be ``full``, ``intra_molecular``,
        or ``inter_molecular``.
    number_i : int, default 0
        Number of atom of the target group to consider for the calculation.
        If ``number_i = 0``, all atoms are considered.
    order : str, default ``m0``
        Order of the analysis, which can be ``m0`` or ``m012``.
        With ``m0``, only the spherical harmonic of order 0 is considered, 
        this is only valid for isotropic systems. For nonisotropic systems,
        use ``m012``.
    f0 : int, default ``None``
        Frequency at which ``T1`` and ``T2`` are calculated.
        If ``None``, ``f = 0`` is used.
    actual_dt : float, default ``None``
        Can be used to specify a different time interface between frames than the 
        one detected by MDAnalysis.
    hydrogen_per_atom : float, default 1.0
        Specify the number of hydrogen per atom, usefull for 
        coarse-grained simulations.
    start : int, default 0
        To specify the first frame to be considered for the analysis.
    stop : int, default 0
        To specify the last frame to be considered for the analysis.
    step : int, default 1
        To jump frame during the analysis.  
    pdb : bool, default True
        To turn off/on the periodic boundary condition treatment.            
    """

    def __init__(self,
                 u: MDAnalysis.Universe,
                 atom_group: MDAnalysis.AtomGroup,
                 type_analysis: str = "full",
                 number_i: int =0,
                 order: str ="m0",
                 f0: float =None,
                 actual_dt: float =None,
                 hydrogen_per_atom: float =1.0,
                 spin: float =1/2,
                 start: int = 0,
                 stop: int = 0,
                 step: int = 1,
                 pbc: bool = True
                 ):
        """Initialise class NMR.
        """

        self.u = u
        if len(atom_group) == 0:
            raise ValueError("Missing atom group")
        elif len(atom_group) == 1:
            self.target_i = atom_group[0]
            self.neighbor_j = atom_group[0]
        elif len(atom_group) == 2:
            self.target_i = atom_group[0]
            self.neighbor_j = atom_group[1]
        elif len(atom_group) >= 3:
            raise ValueError("More than 2 atom groups")
        self.type_analysis = type_analysis
        self.number_i = number_i
        self.order = order
        self.actual_dt = actual_dt
        self.hydrogen_per_atom = hydrogen_per_atom
        self.spin = spin
        self.start = start
        if stop == 0:
            self.stop = u.trajectory.n_frames
        else:
            self.stop = stop
        self.step = step
        self.pbc = pbc
        self.data = None
        self.rij = None
        self.r = None
        self.theta = None
        self.phi = None
        self.sh20 = None
        self.sh21 = None
        self.sh22 = None
        self.gij = None
        self.T1 = None
        self.T2 = None
        self.tau = None
        self.delta_omega = None

        self._verify_entry()
        self._define_constants()
        self._read_universe()
        self._select_target_i()

        # Loop on all the atom of group i
        for _cpt_i, _ in enumerate(self.index_i):
            self._cpt_i = _cpt_i
            #self.atom_index_i = atom_index_i # not used
            self._select_atoms_group_i()
            self._select_atoms_group_j()
            if _cpt_i == 0:
                self._initialise_data()
            self._evaluate_correlation_ij()

        self._calculate_fourier_transform()
        self._calculate_spectrum()
        self.f0 = f0
        self._calculate_relaxationtime()
        self._calculate_tau()
        self._calculate_secondmoment()

    def _verify_entry(self):
        """Verify that entries are correct."""
        _possible_analysis = ["inter_molecular", "intra_molecular", "full"]
        if self.type_analysis not in _possible_analysis:
            raise ValueError("type_analysis can only be inter_molecular, intra_molecular, and full.")
        _possible_order = ["m0", "m012"]
        if self.order not in _possible_order:
            raise ValueError("order can only be m0 and m012.")

    def _define_constants(self):
        """Define the pre-factor K.

        Gamma is the gyromagnetic constant in Hz/T, and
        K has the units of m^6/s^2
        # @tocheck the units
        # @tofix do atoms have all the same gyromag ratio?
        """
        self._GAMMA = 2 * np.pi * 42.6e6
        self._K = (3 / 2) * (cst.mu_0 / 4 / np.pi) ** 2 \
            * cst.hbar ** 2 * self._GAMMA ** 4 * self.spin * (1 + self.spin)

    def _read_universe(self):
        """Create atom groups from chosen atom selections.

        The created groups must not be empty.
        """
        self.group_target_i = self.u.select_atoms(self.target_i)
        self.group_neighbor_j = self.u.select_atoms(self.neighbor_j)
        # Assert that both groups contain atoms
        if self.group_target_i.atoms.n_atoms == 0:
            raise ValueError("Empty atom groups i")
        elif self.group_neighbor_j.atoms.n_atoms == 0:
            raise ValueError("Empty atom groups j")

    def _select_target_i(self):
        if self.number_i == 0:
            self.index_i = np.array(self.group_target_i.atoms.indices)
        elif self.number_i > self.group_target_i.atoms.n_atoms:
            print('Note : number_i is larger than the number of atoms in group target i\n'
                  '-> All the atoms of the group i have been selected')
            self.index_i = np.array(self.group_target_i.atoms.indices)
        else:
            self.index_i = np.array(random.choices(self.group_target_i.atoms.indices, k=self.number_i))
    def _select_atoms_group_i(self):
        """Select atoms of the group i for the calculation."""
        self.group_i = self.u.select_atoms('index ' + str(self.index_i[self._cpt_i]))
        self._resids_i = self.group_i.resids[self.group_i.atoms.indices == self.index_i[self._cpt_i]]

    def _select_atoms_group_j(self):
        """Select atoms of the group j for the calculation.

        For intra molecular analysis, group j are made of atoms of the
        same residue as group i.
        For inter molecular analysis, group j are made of atoms of
        different residues as group i.
        For full analysis, group j are made of atoms that are not in group i.
        """
        if self.type_analysis == "intra_molecular":
            _same_residue : bool = self.group_neighbor_j.resids == self._resids_i
            _different_atom : bool = self.group_neighbor_j.indices != self.index_i[self._cpt_i]
            _index_j = self.group_neighbor_j.atoms.indices[_same_residue & _different_atom]
            _str_j = ' '.join(str(e) for e in _index_j)
        elif self.type_analysis == "inter_molecular":
            _different_residue : bool = self.group_neighbor_j.resids != self._resids_i
            _index_j = self.group_neighbor_j.atoms.indices[_different_residue]
            _str_j = ' '.join(str(e) for e in _index_j)
        elif self.type_analysis == "full":
            _different_atom : bool = self.group_neighbor_j.indices != self.index_i[self._cpt_i]
            _index_j = self.group_neighbor_j.atoms.indices[_different_atom]
            _str_j = ' '.join(str(e) for e in _index_j)
        if len(_str_j) == 0:
            raise ValueError("Empty atom groups j \n"
                             "Wrong combination of type_analysis and group selection?")
        else:
            self.group_j = self.u.select_atoms('index ' + _str_j)

    def _initialise_data(self):
        """Initialise data arrays.

        # @tofix if step>1, matrix will be too large
        # either adapt or cut it before calculating correlation
        """
        if self.order == "m0":
            self._data = np.zeros((1, self.u.trajectory.n_frames,
                                  self.group_j.atoms.n_atoms),
                                  dtype=np.float32)
            self.gij = np.zeros((1,  self.u.trajectory.n_frames),
                                 dtype=np.float32)
        elif self.order == "m012":
            self._data = np.zeros((3, self.u.trajectory.n_frames,
                                  self.group_j.atoms.n_atoms),
                                  dtype=np.complex64)
            self.gij = np.zeros((3, self.u.trajectory.n_frames),
                                dtype=np.float32)
        if self.actual_dt is None:
            _timestep = np.round(self.u.trajectory.dt, 4)
        else:
            _timestep = self.actual_dt
        if self.step > 1:
            _timestep *= self.step
        self.t = np.arange(self.u.trajectory.n_frames) * _timestep

    def _evaluate_correlation_ij(self):
        """Evaluate the correlation function.

        Run over the MDA trajectory. If start, stop, or step are
        specified, only a sub-part of the trajectory is analysed.

        Note: if step>1 is given, the timestep is adapted.
        """
        for _cpt, _ts in enumerate(self.u.trajectory[self.start:self.stop:self.step]):
            self._position_i = self.group_i.atoms.positions
            self._position_j = self.group_j.atoms.positions
            self._box = _ts.dimensions
            # ensure that the box is orthonormal
            if not np.all(self._box[3:] == self._box[3:][0]) & np.all(self._box[3:][0] == 90.0):
                raise ValueError("NMRforMD does not accept non-orthogonal box"
                                 "Use triclinic_to_orthorhombic from the package lipyphilic"
                                 "to convert the trajectory file.")
            self._vector_ij()
            self._cartesian_to_spherical()
            self._spherical_harmonic()
            if self.order == 'm0':
                self._data[:, _cpt] = self._sh20
            elif self.order == 'm012':
                self._data[:, _cpt] = self._sh20, self._sh21, self._sh22
        self._calculate_correlation_ij()

    def _calculate_correlation_ij(self):
        """Calculate the correlation function from ."""
        for _idx_j in range(self.group_j.atoms.n_atoms):
            if self.order == 'm0':
                self.gij += correlation_function(self._data[0, :, _idx_j])
            elif self.order == 'm012':
                for _m in range(3):
                    self.gij[_m] += correlation_function(self._data[_m, :, _idx_j])
        self.gij = np.real(self.gij)

    def _calculate_fourier_transform(self):
        # normalise gij by the number of iteration
        self.gij /= self._cpt_i+1
        # dimensionalize the correlation function
        # from A-6 to s2/m6
        # @tocheck units
        self.gij *= self._K / cst.angstrom ** 6
        # for coarse grained models, possibly more than 1 hydrogen per atom
        self.gij *= np.float32(self.hydrogen_per_atom)
        if self.order == 'm0':
            fij = fourier_transform(np.vstack([self.t, self.gij]).T)
            self.f = np.real(fij.T[0])
            self.J = np.real(fij.T[1])
        elif self.order == 'm012':
            for _m in range(3):
                fij = fourier_transform(np.vstack([self.t, self.gij[_m]]).T)
                self.f = np.real(fij.T[0])
                if _m == 0:
                    _J_0 = np.real(fij.T[1])
                elif _m == 1:
                    _J_1 = np.real(fij.T[1])
                elif _m == 2:
                    _J_2 = np.real(fij.T[1])
            self.J = np.array([_J_0, _J_1, _J_2])


    def _vector_ij(self):
        """Calculate distance between position_i and position_j."""
        if self.pbc:
            self._rij = (np.remainder(self._position_i - self._position_j
                         + self._box[:3]/2., self._box[:3]) - self._box[:3]/2.).T
        else:
            self._rij = (self._position_i - self._position_j).T

    def _cartesian_to_spherical(self):
        """Convert cartesian coordinate to spherical."""
        self._r = np.sqrt(self._rij[0]**2 + self._rij[1]**2 + self._rij[2]**2)
        self._theta = np.arctan2(np.sqrt(self._rij[0]**2 + self._rij[1]**2), self._rij[2])
        self._phi = np.arctan2(self._rij[1], self._rij[0])

    def _spherical_harmonic(self):
        """Evaluate spherical harmonic functions from rij vector.

        convention : theta = polar angle, phi = azimuthal angle
        note: scipy uses the opposite convention
        """
        _a_0 = np.sqrt(16 * np.pi / 5)
        _a_1 = np.sqrt(8 * np.pi / 15)
        _a_2 = np.sqrt(32 * np.pi / 15)
        self._sh20 = _a_0 * sph_harm(0, 2, self._phi, self._theta) / np.power(self._r, 3)
        if self.order == "m0":
            self._sh20 = self._sh20.real
        if self.order == "m012":
            self._sh21 = _a_1 * sph_harm(1, 2, self._phi, self._theta) / np.power(self._r, 3)
            self._sh22 = _a_2 * sph_harm(2, 2, self._phi, self._theta) / np.power(self._r, 3)

    def _calculate_spectrum(self):
        """Calculate spectrums R1 and R2 from J."""
        if self.order == "m0":
            _inter1d = interp1d(self.f, self.J, fill_value="extrapolate")
            self.R1 = (_inter1d(self.f)
                       + 4 * _inter1d(2 * self.f))/6
            self.R2 = (3/2*_inter1d(self.f[0])
                       + 5/2*_inter1d(self.f)
                       + _inter1d(2 * self.f))/6
        elif self.order == "m012":
            _inter1d_0 = interp1d(self.f, self.J[0], fill_value="extrapolate")
            _inter1d_1 = interp1d(self.f, self.J[1], fill_value="extrapolate")
            _inter1d_2 = interp1d(self.f, self.J[2], fill_value="extrapolate")
            self.R1 = _inter1d_1(self.f) + _inter1d_2(2 * self.f)
            self.R2 = (1/4)*(_inter1d_0(self.f[0])+ 10*_inter1d_1(self.f)
                             + _inter1d_2(2 * self.f))

    def _calculate_relaxationtime(self):
        """Calculate the relaxation time at a given frequency f0 (default is 0)"""
        if self.f0 is None:
            self.T1 = 1/self.R1[0]
            self.T2 = 1/self.R2[0]
        else:
            idx = find_nearest(self.f, self.f0)
            self.T1 = 1 / self.R1[idx]
            self.T2 = 1 / self.R2[idx]

    def _calculate_tau(self):
        """
        Calculate correlation time using tau = 0.5 J(0) / G(0).

        The unit are in picosecond. If only the 0th m order is used,
        one value for tau is returned, if all three m orders are used,
        three values for tau are returned.
        """
        if self.order == "m0":
            self.tau = 0.5 * (self.J[0] / self.gij[0][0])
        elif self.order == "m012":
            self.tau = np.array([0.5*(self.J[0][0] / self.gij.T[0][0]),
                                 0.5*(self.J[1][0] / self.gij.T[0][1]),
                                 0.5*(self.J[2][0] / self.gij.T[0][2])])
        self.tau /= cst.pico

    def _calculate_secondmoment(self):
        """
        Calculate second moment Delta omega.

        The unit of Delta omega are kHz /2 pi. The formula is G(0) = (1/3)
        (Delta omega)**2 where G(0) is the correlation function
        at time 0. If only the 0th m order is used, one value for Delta omega
        is returned, if all three m orders are used, three values for Delta
        omega are returned.
        """
        if self.order == "m0":
            self.delta_omega = np.sqrt(3*self.gij[0][0])
        elif self.order == "m012":
            self.delta_omega = np.array([np.sqrt(3*self.gij[0][0]),
                                         np.sqrt(3*self.gij[0][1]),
                                         np.sqrt(3*self.gij[0][2])])
        self.delta_omega /= 1000*2*np.pi
