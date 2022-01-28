import MDAnalysis as mda
import nmrformd as NMR
import numpy as np


def test_nmr():
    u = mda.Universe("../tests/bulk_h2o/topology.tpr", "../tests/bulk_h2o/trajectory.xtc")

    gi = "type HW"
    gj = "type HW"
    t = "full"
    N = 5
    nmr_result = NMR.NMR(u, gi, gj, t, number_i=N)

    print(nmr_result.delta_omega)