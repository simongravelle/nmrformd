import MDAnalysis as mda
import nmrformd as NMR
import numpy as np

print(pwd)

#def test_nmr():
u = mda.Universe("bulk_h2o/topology.tpr", "bulk_h2o/trajectory.xtc")

gi = "type HW"
gj = "type HW"
t = "full"
nmr_result = NMR.NMR(u, gi, gj, t, 1, "m012")

print(nmr_result.delta_omega)


