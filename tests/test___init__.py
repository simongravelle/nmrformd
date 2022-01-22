import MDAnalysis as mda
import nmrformd as NMR
import numpy as np

def test_nmr():
    u = mda.Universe("../tests/bulk_h2o/topology.tpr", "../tests/bulk_h2o/trajectory.xtc")

    gi = "type HW"
    gj = "type HW"
    t = "full"
    nmr_result = NMR.NMR(u, gi, gj, t, 100, "m0")

    print(1/nmr_result.R1[0])
    print(1/nmr_result.R2[0])

    assert np.isclose(1/nmr_result.R1[0], 1/nmr_result.R2[0], rtol=1e-4, atol=0)