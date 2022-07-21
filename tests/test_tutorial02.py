import MDAnalysis as mda
import numpy as np

import nmrformd as nmrmd



def test_nmr():
    """Test NMR module using bulk water trajectory."""
    u = mda.Universe("peg_water/PEG_H2O.data",
                     "peg_water/traj.xtc")

    print()
    print()
    n_water_molecules = u.atoms.select_atoms("type 1").atoms.n_atoms
    print(f"The number of water molecules is {n_water_molecules}")
    n_atom_peg = u.atoms.select_atoms("type 3 4 5 6 7").atoms.n_atoms
    print(f"The number of atoms in the peg molecules is {n_atom_peg}")
    timestep = np.int32(u.trajectory.dt)
    print(f"The timestep is {timestep} ps")
    total_time = np.int32(u.trajectory.totaltime)
    print(f"The total simulation time is {total_time} ps")
    print()
    print()
    print()
    print()

    group_H2O = "type 2"
    group_PEG = "type 5 7"
    group_ALL = "type 2 5 7"

    print(u.select_atoms(group_PEG).atoms.indices)
    print(u.select_atoms(group_PEG).atoms.resids)

    nmr_PEG = nmrmd.NMR(u, [group_PEG, group_PEG], number_i=40, type_analysis="intra_molecular")
    print(nmr_PEG.T1)


    print()
    print()
    print()
    print()