import MDAnalysis as mda
import nmrformd as NMR

def test_nmr():

    u = mda.Universe("../test/run.tpr", "../test/run.xtc")

    nmr_result = NMR.NMR(u, "type HW","type HW",1,1)
