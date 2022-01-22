import MDAnalysis as mda
import nmrformd as NMR
import numpy as np
import matplotlib.pyplot as plt

def test_nmr():
    u = mda.Universe("../test/bulk_h2o/topology.tpr", "../test/bulk_h2o/trajectory.xtc")

    nmr_result = NMR.NMR(u, "type HW", "type HW", "full", 1000, "m012")

    print()
    print()
    # print("result "+str(nmr_result.index_i))
    # print("result " + str(nmr_result.group_i))
    # print(nmr_result.group_i.atoms.types[0])
    # print(nmr_result.target_i)
    #print(nmr_result.gij.t)
    #print(nmr_result.gij.t)
    #print(nmr_result.f.shape)
    print(1/nmr_result.R1[0])
    print(1/nmr_result.R2[0])

    #np.savetxt('test.dat', np.stack(nmr_result.t.T,nmr_result.gij[0].T), delimiter=' ')


    print()



