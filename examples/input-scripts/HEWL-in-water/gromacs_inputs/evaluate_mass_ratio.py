import MDAnalysis as mda
import numpy as np

u = mda.Universe("conf.gro")
# Water
water = u.select_atoms("name OW HW1 HW2")
mass_water = np.sum(water.masses)

# Protein
all_name = ' '
for name in np.unique(u.atoms.names):
    if (name != 'OW') & (name != 'HW2') & (name != 'HW1'):
        all_name += name + ' '
protein = u.select_atoms('name '+all_name)
mass_protein = np.sum(protein.masses)

print("mass_water/mass_protein = ", np.round(mass_water/mass_protein,3))