import MDAnalysis as mda
import nmrformd as NMR
import numpy as np
import matplotlib.pyplot as plt

def test_nmr():
    u = mda.Universe("bulk_h2o/topology.tpr", "bulk_h2o/trajectory.xtc")

    nmr_result = NMR.NMR(u, "type HW", "type HW", "full", 100, "m0")

    print(1/nmr_result.R1[0])
    print(1/nmr_result.R2[0])


