import MDAnalysis as mda
import nmrformd as NMR


def test_nmr():
    u = mda.Universe("../test/bulk_h2o/topology.tpr", "../test/bulk_h2o/trajectory.xtc")

    nmr_result = NMR.NMR(u, "type HW", "type HW", "full", 1)

    print()
    print()
    # print("result "+str(nmr_result.index_i))
    # print("result " + str(nmr_result.group_i))
    # print(nmr_result.group_i.atoms.types[0])
    # print(nmr_result.target_i)
    print(nmr_result.group_i)
    print(nmr_result.group_j)
    print()


print()
