import MDAnalysis as mda
import nmrformd as NMR
import numpy as np


def test_nmr():
    u = mda.Universe("../tests/bulk_h2o/topology.tpr", "../tests/bulk_h2o/trajectory.xtc")

    gi = "type HW"
    gj = "type HW"
    t = "full"
    nmr_result = NMR.NMR(u, gi, gj, t, 1, "m0")

    print(nmr_result.calculateT1T2(f0=0))
    print(nmr_result.T1)

    print(nmr_result.f[:100])
    # assert np.isclose(1/nmr_result.R1[0], 1/nmr_result.R2[0], rtol=1e-4, atol=0)