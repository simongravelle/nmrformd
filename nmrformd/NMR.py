 #!/usr/bin/env python3

from scipy import constants as cst
from scipy.special import sph_harm
from scipy.interpolate import interp1d
import numpy as np
import random

from .utilities import fourier_transform, correlation_function, find_nearest


class NMR:
    """Calculate NMR relaxation time from MDAnalysis universe

    Parameters
    ----------

    u : MDAnalysis.core.universe.Universe
        MDAnalysis universe containing all the information describing
        the molecular dynamics system.
    target_i : AtomGroup
        Target group of the atoms ``i`` for the NMR calculation.
    neighbor_j : AtomGroup
        Neighbor group of the atoms ``j`` for the NMR calculation.
    type_analysis : str
        Type of analysis, which can be ``full``, ``intra_molecular``, or ``inter_molecular``.
    number_i : int, default 0
        Number of atom of the group ``target_i`` to consider. If ``number = 0``, all atoms
        are considered.
    order : str, default ``m0``
        Order of the analysis, which can be ``m0`` or ``m012``.
    f0 : int, default ``None``
        Frequency for the calculation of ``T1`` and ``T2``.  If ``None``, ``f = 0`` is used.
    """

    def __init__(self,
                 u,
                 target_i,
                 neighbor_j,
                 type_analysis,
                 number_i = 0,
                 order = "m0",
                 f0 = None,
                 actual_dt = 0,
                 hydrogen_per_atom = 1.0,
                 spin = 1/2
                 ):

        self.u = u
        self.target_i = target_i
        self.neighbor_j = neighbor_j
        self.type_analysis = type_analysis
        self.number_i = number_i
        self.order = order
        self.actual_dt = actual_dt
        self.hydrogen_per_atom = hydrogen_per_atom
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
        self.spin = spin

        self._define_constants()
        self._read_universe()

        self._select_target_i()

        for cpt_i, i in enumerate(self.index_i):
            self.cpt_i = cpt_i
            self.i = i
            self._select_proton()
            if cpt_i == 0:
                self._initialise_data()
            self._evaluate_correlation_ij()

        self._calculate_fourier_transform()
        self._calculate_spectrum()
        self.f0 = f0
        self._calculate_relaxationtime(self.f0)
        self._calculate_tau()
        self._calculate_secondmoment()

    def _define_constants(self):
        self.GAMMA = 2*np.pi*42.6e6  # gyromagnetic constant in Hz/T
        #self.K = (3*np.pi/5)*(cst.mu_0/4/np.pi)**2*cst.hbar**2*self.GAMMA**4  # m6/s2
        self.K = (3 / 2) * (cst.mu_0 / 4 / np.pi) ** 2 \
                 * cst.hbar ** 2 * self.GAMMA ** 4 * self.spin * (1 + self.spin)  # m6/s2

    def _read_universe(self):
        self.group_target_i = self.u.select_atoms(self.target_i)
        assert self.group_target_i.atoms.n_atoms > 0, "empty target group i"
        self.group_neighbor_j = self.u.select_atoms(self.neighbor_j)
        assert self.group_neighbor_j.atoms.n_atoms > 0, "empty neighbor group j"

    def _select_target_i(self):
        if self.number_i == 0:
            self.index_i = np.array(self.group_target_i.atoms.indices)
        elif self.number_i > self.group_target_i.atoms.n_atoms:
            print('number_i larger than the number of atom in group target i, all atoms selected')
            self.index_i = np.array(self.group_target_i.atoms.indices)
        else:
            self.index_i = np.array(random.choices(self.group_target_i.atoms.indices, k=self.number_i))

    def _select_proton(self):
        assert (self.type_analysis == "inter_molecular") | \
               (self.type_analysis == "intra_molecular") | \
               (self.type_analysis == "full"), \
               'unknown value for type_analysis, choose inter_molecular or intra_molecular or full'

        self.group_i = self.u.select_atoms('index '+str(self.index_i[self.cpt_i]))
        self.resids_i = self.group_i.resids[self.group_i.atoms.indices == self.index_i[self.cpt_i]]

        if self.type_analysis == "inter_molecular":
            self.index_j = \
                self.group_neighbor_j.atoms.indices[(self.group_neighbor_j.resids == self.resids_i)
                                                    & (self.group_neighbor_j.indices != self.index_i[self.cpt_i])]
            self.str_j = ' '.join(str(e) for e in self.index_j)
            self.group_j = self.u.select_atoms('index ' + self.str_j)
        elif self.type_analysis == "intra_molecular":
            self.index_j = self.group_neighbor_j.atoms.indices[self.group_neighbor_j.resids != self.resids_i]
            self.str_j = ' '.join(str(e) for e in self.index_j)
            self.group_j = self.u.select_atoms('index ' + self.str_j)
        elif self.type_analysis == "full":
            self.index_j = \
                self.group_neighbor_j.atoms.indices[self.group_neighbor_j.indices != self.index_i[self.cpt_i]]
            self.str_j = ' '.join(str(e) for e in self.index_j)
            self.group_j = self.u.select_atoms('index ' + self.str_j)

    def _initialise_data(self):
        assert (self.order == "m0") | (self.order == "m012"), 'unknown value for  order, choose m0 or m012'
        if self.order == "m0":
            self.data = np.zeros((1, self.u.trajectory.n_frames, self.group_j.atoms.n_atoms), dtype=np.float32)
            self.gij = np.zeros((1,  self.u.trajectory.n_frames), dtype=np.float32)
        elif self.order == "m012":
            self.data = np.zeros((3, self.u.trajectory.n_frames, self.group_j.atoms.n_atoms), dtype=np.complex64)
            self.gij = np.zeros((3, self.u.trajectory.n_frames), dtype=np.float32)
        if self.actual_dt == 0:
            self.t = np.arange(self.u.trajectory.n_frames) * np.round(self.u.trajectory.dt, 4)
        else :
            self.t = np.arange(self.u.trajectory.n_frames) * np.round(self.actual_dt, 4)

    def _evaluate_correlation_ij(self):
        for cpt, ts in enumerate(self.u.trajectory):
            self.position_i = self.group_i.atoms.positions
            self.position_j = self.group_j.atoms.positions
            self.box = self.u.dimensions
            self._vector_ij()
            self._cartesian_to_spherical()
            self._spherical_harmonic()
            if self.order == 'm0':
                self.data[:, cpt] = self.sh20
            elif self.order == 'm012':
                self.data[:, cpt] = self.sh20, self.sh21, self.sh22
        self._calculate_correlation_ij()

    def _calculate_correlation_ij(self):
        for idx_j in range(self.group_j.atoms.n_atoms):
            if self.order == 'm0':
                self.gij += correlation_function(self.data[0, :, idx_j])
            elif self.order == 'm012':
                for m in range(3):
                    self.gij[m] += correlation_function(self.data[m, :, idx_j])
        self.gij = np.real(self.gij)

    def _calculate_fourier_transform(self):
        self.gij /= self.cpt_i+1
        self.gij *= self.K / cst.angstrom ** 6
        self.gij *= np.float32(self.hydrogen_per_atom)
        if self.order == 'm0':
            fij = fourier_transform(np.vstack([self.t, self.gij]).T)
            self.f = np.real(fij.T[0])
            self.J_0 = np.real(fij.T[1])

        elif self.order == 'm012':
            for m in range(3):
                fij = fourier_transform(np.vstack([self.t, self.gij[m]]).T)
                self.f = np.real(fij.T[0])
                if m == 0:
                    self.J_0 = np.real(fij.T[1])
                elif m == 1:
                    self.J_1 = np.real(fij.T[1])
                elif m == 2:
                    self.J_2 = np.real(fij.T[1])

    def _vector_ij(self):
        """calculate 3D vector between position_i and position_j assuming pbc"""
        self.rij = (np.remainder(self.position_i - self.position_j + self.box[:3]/2., self.box[:3]) - self.box[:3]/2.).T

    def _cartesian_to_spherical(self):
        self.r = np.sqrt(self.rij[0]**2 + self.rij[1]**2 + self.rij[2]**2)
        self.theta = np.arctan2(np.sqrt(self.rij[0]**2 + self.rij[1]**2), self.rij[2])
        self.phi = np.arctan2(self.rij[1], self.rij[0])

    def _spherical_harmonic(self):
        """evaluate harmonic functions from rij vector
        convention : theta = polar angle, phi = azimuthal angle
        note: scipy uses the opposite convention
        """
        self.sh20 = np.sqrt(16 * np.pi/ 5) * sph_harm(0, 2, self.phi, self.theta) / np.power(self.r, 3)
        if self.order == 'm012':
            self.sh21 = np.sqrt(8 * np.pi/ 15) * sph_harm(1, 2, self.phi, self.theta) / np.power(self.r, 3)
            self.sh22 = np.sqrt(32 * np.pi/ 15) * sph_harm(2, 2, self.phi, self.theta) / np.power(self.r, 3)

    def _calculate_spectrum(self):
        interpolation_0 = interp1d(self.f, self.J_0, fill_value="extrapolate")
        if self.order == "m0":
            self.R1 = (interpolation_0(self.f) + 4 * interpolation_0(2 * self.f))/6
            self.R2 = (3/2*interpolation_0(self.f[0])+5/2*interpolation_0(self.f) + interpolation_0(2 * self.f))/6
        elif self.order == "m012":
            interpolation_1 = interp1d(self.f, self.J_1, fill_value="extrapolate")
            interpolation_2 = interp1d(self.f, self.J_2, fill_value="extrapolate")
            self.R1 = interpolation_1(self.f) + interpolation_2(2 * self.f)
            self.R2 = (1/4)*(interpolation_0(self.f[0])+10*interpolation_1(self.f) + interpolation_2(2 * self.f))

    def _calculate_relaxationtime(self, f0=None):
        if f0 is None:
            self.T1 = 1/self.R1[0]
            self.T2 = 1/self.R2[0]
        else:
            idx = find_nearest(self.f, f0)
            self.T1 = 1 / self.R1[idx]
            self.T2 = 1 / self.R2[idx]

    def _calculate_tau(self):
        """
        Calculate correlation time using tau = 0.5 J(0) / G(0).

        The unit are in picosecond. If only the 0th m order is used, one value for tau
        is returned, if all three m orders are used, three values for tau are used.
        """
        if self.order == "m0":
            self.tau = 0.5 * (self.J_0[0] / self.gij[0][0])
        elif self.order == "m012":
            self.tau = np.array([0.5*(self.J_0[0] / self.gij[0][0]),
                                 0.5*(self.J_1[0] / self.gij[0][1]),
                                 0.5*(self.J_2[0] / self.gij[0][2])])
        self.tau /= cst.pico

    def _calculate_secondmoment(self):
        """
        Calculate second moment Delta omega.

        The unit of Delta omega are kHz /2 pi. The formula is G(0) = (1/3)
        (Delta omega)**2 where G(0) is the correlation function
        at time 0. If only the 0th m order is used, one value for Delta omega
        is returned, if all three m orders are used, three values for Delta omega are used.
        """
        if self.order == "m0":
            self.delta_omega = np.sqrt(3*self.gij[0][0])
        elif self.order == "m012":
            self.delta_omega = np.array([np.sqrt(3*self.gij[0][0]),
                                         np.sqrt(3*self.gij[0][1]),
                                         np.sqrt(3*self.gij[0][2])])
        self.delta_omega /= 1000*2*np.pi
