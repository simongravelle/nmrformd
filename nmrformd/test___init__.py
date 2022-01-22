import MDAnalysis as mda
import nmrformd as NMR


def test_nmr():
    u = mda.universe("../tests/bulk_h2o/topology.tpr", "../tests/bulk_h2o/trajectory.xtc")

    nmr_result = NMR.NMR(u, "type HW", "type HW", "full", 100, "m0")

    print(1/nmr_result.R1[0])
    print(1/nmr_result.R2[0])




