
import numpy as np

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
        self._read_universe()

    def _read_universe(self):
        self.group_target_i = self.u.select_atoms(self.target_i)
        self.group_neighbor_j = self.u.select_atoms(self.neighbor_j)

        print(self.group_target_i.atoms.n_atoms)
        print(self.group_neighbor_j.atoms.n_atoms)




