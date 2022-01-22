from scipy import constants as cst
import numpy as np
import random


class NMR:
    def __init__(self,
                 u,
                 target_i,
                 neighbor_j,
                 type_analysis,
                 number_j):
        self.u = u
        self.target_i = target_i
        self.neighbor_j = neighbor_j
        self.type_analysis = type_analysis
        self.number_j = number_j
        self._define_constants()
        self._read_universe()

    def _define_constants(self):
        self.gamma = 2*np.pi*42.6e6 # gyromagnetic constant in Hz/T
        self.K = (3*np.pi/5)*(cst.mu_0/4/np.pi)**2*(cst.hbar)**2*(self.gamma)**4 # m6/s2

    def _read_universe(self):
        self.group_target_i = self.u.select_atoms(self.target_i)
        assert self.group_target_i.atoms.n_atoms > 0, "empty target group i"
        self.group_neighbor_j = self.u.select_atoms(self.neighbor_j)
        assert self.group_neighbor_j.atoms.n_atoms > 0, "empty neighbor group j"
        self._select_proton()

    def _select_proton(self):
        self.index_i = np.array(random.choices(self.group_target_i.atoms.indices, k=1))
        self.group_i = self.u.select_atoms('index '+str(self.index_i[0]))
        self.resids_i = self.group_i.resids[self.group_i.atoms.indices == self.index_i[0]]

        if self.type_analysis == "inter_molecular":
            self.index_j = self.group_neighbor_j.atoms.indices[(self.group_neighbor_j.resids == self.resids_i)
                                                               & (self.group_neighbor_j.indices != self.index_i[0])]
            self.str_j = ' '.join(str(e) for e in self.index_j)
            self.group_j = self.u.select_atoms('index ' + self.str_j)
        elif self.type_analysis == "intra_molecular":
            self.index_j = self.group_neighbor_j.atoms.indices[self.group_neighbor_j.resids != self.resids_i]
            self.str_j = ' '.join(str(e) for e in self.index_j)
            self.group_j = self.u.select_atoms('index ' + self.str_j)
        elif self.type_analysis == "full":
            self.index_j = self.group_neighbor_j.atoms.indices[self.group_neighbor_j.indices != self.index_i[0]]
            self.str_j = ' '.join(str(e) for e in self.index_j)
            self.group_j = self.u.select_atoms('index ' + self.str_j)

        self._evaluate_correlation_ij()

    def _evaluate_correlation_ij(self):
        for ts in self.u.trajectory:
            self.position_i = self.group_i.atoms.positions
            self.position_j = self.group_j.atoms.positions
            self.box = self.u.dimensions
