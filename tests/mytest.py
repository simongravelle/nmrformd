import MDAnalysis as mda
import nmrformd as NMR
import numpy as np

def test_nmr():
    u = mda.Universe("bulk_h2o/topology.tpr", "bulk_h2o/trajectory.xtc")

    group_i = "type HW and index 0:20"
    group_j = "type HW and index 0:20"
    analysis_type = "full"
    number_i = 0
    order = "m0"
    nmr_result = NMR.NMR(u, group_i,
                         group_j, analysis_type,
                         number_i, order)
    T1 = nmr_result.T1
    T2 = nmr_result.T2
    Dw = nmr_result.delta_omega
    tau = nmr_result.tau

    assert np.isclose(T1, 1.029287074574358, rtol=1e-4, atol=0)
    assert np.isclose(T2, 1.029287074574358, rtol=1e-4, atol=0)
    assert np.isclose(Dw, 56.996441578600304, rtol=1e-4, atol=0)
    assert np.isclose(tau, 13.635804446316115, rtol=1e-4, atol=0)

    group_i = "type HW and index 0:20"
    group_j = "type HW and index 0:20"
    analysis_type = "inter_molecular"
    number_i = 0
    order = "m0"
    nmr_result = NMR.NMR(u, group_i,
                         group_j, analysis_type,
                         number_i, order)
    T1 = nmr_result.T1
    T2 = nmr_result.T2
    Dw = nmr_result.delta_omega
    tau = nmr_result.tau

    assert np.isclose(T1, 1.0305, rtol=1e-4, atol=0)
    assert np.isclose(T2, 1.0305, rtol=1e-4, atol=0)
    assert np.isclose(Dw, 56.7908, rtol=1e-4, atol=0)
    assert np.isclose(tau, 13.7180, rtol=1e-4, atol=0)

    group_i = "type HW and index 0:40"
    group_j = "type HW and index 0:40"
    analysis_type = "intra_molecular"
    number_i = 0
    order = "m0"
    nmr_result = NMR.NMR(u, group_i,
                         group_j, analysis_type,
                         number_i, order)
    T1 = nmr_result.T1
    T2 = nmr_result.T2
    Dw = nmr_result.delta_omega
    tau = nmr_result.tau

    assert np.isclose(T1, 400.4619, rtol=1e-4, atol=0)
    assert np.isclose(T2, 400.4619, rtol=1e-4, atol=0)
    assert np.isclose(Dw, 5.99748, rtol=1e-4, atol=0)
    assert np.isclose(tau, 3.16528, rtol=1e-4, atol=0)

    group_i = "type HW and index 0:40"
    group_j = "type HW and index 0:40"
    analysis_type = "full"
    number_i = 0
    order = "m012"
    f0 = 10000
    nmr_result = NMR.NMR(u, group_i,
                         group_j, analysis_type, 
                         number_i, order, f0)
    T1 = nmr_result.T1
    T2 = nmr_result.T2
    Dw = nmr_result.delta_omega
    tau = nmr_result.tau

    assert np.isclose(T1, 4.931728, rtol=1e-4, atol=0)
    assert np.isclose(T2, 3.853544, rtol=1e-4, atol=0)
    assert np.isclose(Dw[0], 56.152647, rtol=1e-4, atol=0)
    assert np.isclose(tau[0], 6.891118, rtol=1e-4, atol=0)
