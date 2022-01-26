from scipy import constants as cst
from scipy.special import sph_harm
from scipy.interpolate import interp1d
import numpy as np
import random


def fourier_transform(a):
    """
    Wrapper function that takes the data in real space with
    columns time, frequency and returns the data in Fourier space
    with columns frequency, signal
    Units in : ps, signal
    Units out : MHz, s*signal

    Credit to the fourier transform function of MAICoS
    https://maicos-devel.gitlab.io/maicos/index.html
    """
    dt = (a[1, 0] - a[0, 0]) * cst.pico  # second
    return np.vstack((np.fft.rfftfreq(len(a), dt) / cst.mega, np.fft.rfft(a[:, 1], norm=None) * dt * 2)).T


def correlation_function(a, b=None):
    """Uses fast fourier transforms to calculate the correlation function.

    If two arrays are given, the cross-correlation is returned.
    If one array is given, the auto-correlation is returned.

    Credit to the correlation function of MAICoS
    https://maicos-devel.gitlab.io/maicos/index.html
    """
    if len(a.shape) > 1:
        a2 = np.append(a,
                       np.zeros((2 ** int(np.ceil((np.log(len(a)) / np.log(2)))) - len(a), a.shape[1])),
                       axis=0)
    else:
        a2 = np.append(a, np.zeros(2 ** int(np.ceil((np.log(len(a)) / np.log(2)))) - len(a)), axis=0)
    data_a = np.append(a2, np.zeros(a2.shape), axis=0)
    fra = np.fft.fft(data_a, axis=0)
    if b is None:
        sf = np.conj(fra) * fra
    else:
        if len(b.shape) > 1:
            b2 = np.append(b, np.zeros(
                (2 ** int(np.ceil((np.log(len(b)) / np.log(2)))) - len(b), b.shape[1])), axis=0)
        else:
            b2 = np.append(b, np.zeros(2 ** int(np.ceil((np.log(len(b)) / np.log(2)))) - len(b)),
                           axis=0)
        data_b = np.append(b2, np.zeros(b2.shape), axis=0)
        frb = np.fft.fft(data_b, axis=0)
        sf = np.conj(fra) * frb
    res = np.fft.ifft(sf, axis=0)
    if len(a.shape) > 1:
        # average over all particles/molecules
        cor = (np.real(res[:len(a)]) / np.array(range(len(a), 0, -1))[:, np.newaxis]).mean(axis=1)
    else:
        cor = np.real(res[:len(a)]) / np.array(range(len(a), 0, -1))
    return cor


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


class NMR:
    def __init__(self,
                 u,
                 target_i,
                 neighbor_j,
                 type_analysis,
                 number_i,
                 order):
        self.u = u
        self.target_i = target_i
        self.neighbor_j = neighbor_j
        self.type_analysis = type_analysis
        self.number_i = number_i
        self.order = order
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
        self.calculate_relaxationtime()
        self.calculate_tau()
        self.calculate_secondmoment()

    def _define_constants(self):
        self.gamma = 2*np.pi*42.6e6  # gyromagnetic constant in Hz/T
        self.K = (3*np.pi/5)*(cst.mu_0/4/np.pi)**2*cst.hbar**2*self.gamma**4  # m6/s2

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
        self.t = np.arange(self.u.trajectory.n_frames) * np.round(self.u.trajectory.dt, 4)

    def _evaluate_correlation_ij(self):
        for cpt, ts in enumerate(self.u.trajectory):
            self.position_i = self.group_i.atoms.positions
            self.position_j = self.group_j.atoms.positions
            self.box = self.u.dimensions
            self.vector_ij()
            self.cartesian_to_spherical()
            self.spherical_harmonic()
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

    def vector_ij(self):
        """calculate 3D vector between position_i and position_j assuming pbc"""
        self.rij = (np.remainder(self.position_i - self.position_j + self.box[:3]/2., self.box[:3]) - self.box[:3]/2.).T

    def cartesian_to_spherical(self):
        self.r = np.sqrt(self.rij[0]**2 + self.rij[1]**2 + self.rij[2]**2)
        self.theta = np.arctan2(np.sqrt(self.rij[0]**2 + self.rij[1]**2), self.rij[2])
        self.phi = np.arctan2(self.rij[1], self.rij[0])

    def spherical_harmonic(self):
        """evaluate harmonic functions from rij vector
        convention : theta = polar angle, phi = azimuthal angle
        note: scipy uses the opposite convention
        """
        self.sh20 = sph_harm(0, 2, self.phi, self.theta) / np.power(self.r, 3)
        if self.order == 'm012':
            self.sh21 = sph_harm(1, 2, self.phi, self.theta) / np.power(self.r, 3)
            self.sh22 = sph_harm(2, 2, self.phi, self.theta) / np.power(self.r, 3)

    def _calculate_spectrum(self):
        if self.order == "m0":
            interpolation_0 = interp1d(self.f, self.J_0, fill_value="extrapolate")
            self.R1 = interpolation_0(self.f)+4*interpolation_0(2*self.f)
            self.R2 = 3/2*interpolation_0(self.f[0])+5/2*interpolation_0(self.f) + interpolation_0(2 * self.f)
        elif self.order == "m012":
            interpolation_0 = interp1d(self.f, self.J_0, fill_value="extrapolate")
            interpolation_1 = interp1d(self.f, self.J_1, fill_value="extrapolate")
            interpolation_2 = interp1d(self.f, self.J_2, fill_value="extrapolate")
            self.R1 = interpolation_1(self.f) + interpolation_2(2 * self.f)
            self.R2 = (1/4)*(interpolation_0(self.f[0])+10*interpolation_1(self.f) + interpolation_2(2 * self.f))

    def calculate_relaxationtime(self, f0=None):
        if f0 is None:
            self.T1 = 1/self.R1[0]
            self.T2 = 1/self.R2[0]
        else:
            idx = find_nearest(self.f, f0)
            self.T1 = 1 / self.R1[idx]
            self.T2 = 1 / self.R2[idx]

    def calculate_tau(self):
        """
        Calculate correlation time using tau = 0.5 J(0) / G(0) and convert in pico second
        """
        if self.order == "m0":
            self.tau = 0.5 * (self.J_0[0] / self.gij[0][0])
        elif self.order == "m012":
            self.tau = np.array([0.5*(self.J_0[0] / self.gij[0][0]),
                                 0.5*(self.J_1[0] / self.gij[0][1]),
                                 0.5*(self.J_2[0] / self.gij[0][2])])
        self.tau /= cst.pico

    def calculate_secondmoment(self):
        """
        Calculate second moment Delta omega from G(0) = (1/3) (Delta omega)**2
        and convert in the units of kHz /2 pi
        """
        if self.order == "m0":
            self.delta_omega = np.sqrt(3*self.gij[0][0])
        elif self.order == "m012":
            self.delta_omega = np.array([np.sqrt(3*self.gij[0][0]),
                                         np.sqrt(3*self.gij[0][1]),
                                         np.sqrt(3*self.gij[0][2])])
        self.delta_omega /= 1000*2*np.pi
