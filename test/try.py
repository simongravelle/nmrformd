import MDAnalysis as mda
u = mda.Universe('run.tpr', 'run.xtc')
print(u.atoms.n_atoms)
from nmrformd import NMR

cal_nmr = NMR(u,'HW1','HW2',5)
cal_nmr