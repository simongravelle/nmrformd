import MDAnalysis as mda
import numpy as np

import nmrformd as nmrmd

u = mda.Universe("bulk_h2o/topology.tpr", "bulk_h2o/trajectory.xtc")

group_i = "type HW and index 0:20"
group_j = "type HW and index 0:20"
analysis_type = "full"
number_i = 0
order = "m0"
nmr_result = nmrmd.NMR(u, group_i,
                     group_j, analysis_type,
                     number_i, order)
T1 = nmr_result.T1
T2 = nmr_result.T2
Dw = nmr_result.delta_omega
tau = nmr_result.tau

print()
print(T1)

assert np.isclose(T1, 1.029287074574358, rtol=1e-4, atol=0)
assert np.isclose(T2, 1.029287074574358, rtol=1e-4, atol=0)
assert np.isclose(Dw, 56.996441578600304, rtol=1e-4, atol=0)
assert np.isclose(tau, 13.635804446316115, rtol=1e-4, atol=0)