import sys
import readlammps
import numpy as np
from copy import deepcopy

sys.path.append("./")
from utilities import generate_random_location, generate_random_orientation, \
                  search_closest_neighbor

cpt_atom = 0
cpt_na = 0
cpt_cl = 0
cpt_h2o = 0
atom_XYZ = []
atom_type = []
atom_id = []
atom_charge = []
atom_mass = []
atom_name = []
residue_number = []
residue_name = []

## import silica and counter ions
silica_info = readlammps.ReadLammpsData("silica_surface/silica_Q3_4_7OH_18pct_ion.lammps05", 
                                        verbose=False)
silica_width = np.max(silica_info.coords.T[2])-np.min(silica_info.coords.T[2])

## initial box dimension in Angstrom
Lx = silica_info.box_dimensions[0][1]
Ly = silica_info.box_dimensions[1][1]
Lz = 50

## place silica
name_SiOH = ['O', 'O', 'Si', 'H']
type_SiOH = ['OC23', 'OC24', 'SC4', 'HOY']
mass_SiOH = [15.999, 15.999, 28.086, 1.008]
for cpt, id_atom in enumerate(silica_info.atom_labels):
    if id_atom <= 4: # exclude Na ions
        atom_XYZ.append(silica_info.coords[cpt]+[0, 0, -silica_width/2])
        atom_type.append(type_SiOH[id_atom-1])
        atom_id.append(cpt_atom+1)
        atom_charge.append(silica_info.atom_charges[cpt])
        atom_mass.append(mass_SiOH[id_atom-1])
        atom_name.append(name_SiOH[id_atom-1])
        residue_number.append(1)
        residue_name.append("SiOH")
        cpt_atom += 1

## place counter ions
name_Na = ['Na']
type_Na = ['Na']
mass_Na = [22.99]
for cpt, id_atom in enumerate(silica_info.atom_labels):
    if id_atom == 5: # only Na ions
        atom_XYZ.append(silica_info.coords[cpt]+[0, 0, -silica_width/2])
        atom_type.append('Na')
        atom_id.append(cpt_atom+1)
        atom_name.append('Na')
        residue_number.append(2+cpt_na)
        residue_name.append("Na")
        cpt_atom += 1
        cpt_na += 1

## define water
H2O_XYZ = np.array([[0, 0, 0], \
       [0.5858,   0.757, 0.0], \
       [0.5858,   -0.757,  0.0], \
       [0.104,  0.0, 0.0]])
H2O_type = ['OW', 'HW1', 'HW2', 'MW']

## place water
cut_off = 2 # Angstrom
delta = 2.9
# ordered placement
for x in np.arange(0, Lx, delta):
    for y in np.arange(0, Ly, delta):
        for z in np.arange(0, Lz, delta):
            dx, dy, dz = generate_random_location(0.1, 0.1, 0.1, recenter = True)
            H2O_XYZ = generate_random_orientation(H2O_XYZ)
            min_distance = search_closest_neighbor(np.array(atom_XYZ),
                                                   H2O_XYZ + [x, y, z] + [dx, dy, dz],
                                                   Lx, Ly, Lz)
            if min_distance > cut_off:
                cpt_h2o += 1
                for cpt in range(len(H2O_type)):
                    atom_XYZ.append(H2O_XYZ[cpt]+ [x, y, z])
                    atom_type.append(H2O_type[cpt])
                    atom_id.append(cpt_atom+1)
                    atom_name.append(H2O_type[cpt])
                    residue_number.append(1+cpt_na+cpt_h2o)
                    residue_name.append("SOL")
                    cpt_atom += 1  

print(str(cpt_h2o) + " inserted water molecules")
vol_water = cpt_h2o/6.022e23*0.018 # kg or litter
h_estimated = vol_water/1e3 / (Lx*1e-10)/(Ly*1e-10)*1e9
print("final wall distance ~ " + str(np.round(h_estimated,2))+" nm")
quantity_ion = (cpt_cl+cpt_na)/6.022e23 # mol
concentration_ion = quantity_ion/vol_water # mol/kg
print('salt concentration = '+str(np.round(concentration_ion,2))+' M')

## write conf.gro
f = open('conf.gro', 'w')
f.write('Silica slit with water and counter ions\n')
f.write(str(cpt_atom)+"\n")
for cpt in range(cpt_atom):
    # residue number (5 positions, integer) 
    f.write("{: >5}".format(str(residue_number[cpt])))
    # residue name (5 characters) 
    f.write("{: >5}".format(str(residue_name[cpt])))
    # atom name (5 characters) 
    f.write("{: >5}".format(str(atom_name[cpt])))
    # atom number (5 positions, integer)
    f.write("{: >5}".format(str(atom_id[cpt]))) 
    # position (in nm, x y z in 3 columns, each 8 positions 
    #with 3 decimal places)
    f.write("{: >8}".format(str("{:.3f}".format(atom_XYZ[cpt][0]/10))))
    f.write("{: >8}".format(str("{:.3f}".format(atom_XYZ[cpt][1]/10))))
    f.write("{: >8}".format(str("{:.3f}".format(atom_XYZ[cpt][2]/10))))
    f.write("\n")
f.write("{: >10}".format(str("{:.5f}".format(Lx/10))))
f.write("{: >10}".format(str("{:.5f}".format(Ly/10))))
f.write("{: >10}".format(str("{:.5f}".format(Lz/10))))
f.write("\n")
f.close()

## write topol.top
f = open('topol.top', 'w')
f.write('#include "ff/forcefield.itp"\n')
f.write('#include "ff/silica.itp"\n')
f.write('#include "ff/tip4peps.itp"\n')
f.write('#include "ff/ions.itp"\n\n')
f.write('[ System ]\n')
f.write('Silica slit with water and counter ions\n\n')
f.write('[ Molecules ]\n\n')
f.write('SiOH 1\n')
if cpt_na>0:
    f.write('NA '+ str(cpt_na)+'\n')
if cpt_cl>0:
    f.write('CL '+ str(cpt_cl)+'\n')
f.write('SOL '+ str(cpt_h2o)+'\n')
f.close()

## re-number bond and angle
prev_bond = silica_info.bonds
new_bond = deepcopy(prev_bond)
new_id = 0
for cpt, old_id in enumerate(silica_info.atom_ids):
    label = silica_info.atom_labels[cpt]
    if label <= 4:
        new_id += 1
        where_to_replace = np.where(prev_bond.T[0] == old_id)
        if len(where_to_replace[0]) > 0:
            for line in where_to_replace[0]:
                new_bond[line,0] = new_id
        where_to_replace = np.where(prev_bond.T[1] == old_id)
        if len(where_to_replace[0]) > 0:
            for line in where_to_replace[0]:
                new_bond[line,1] = new_id 
silica_info.bonds = new_bond

prev_angle = silica_info.angles
new_angle = deepcopy(prev_angle)
new_id = 0
for cpt, old_id in enumerate(silica_info.atom_ids):
    label = silica_info.atom_labels[cpt]
    if label <= 4:
        new_id += 1
        where_to_replace = np.where(prev_angle.T[0] == old_id)
        if len(where_to_replace[0]) > 0:
            for line in where_to_replace[0]:
                new_angle[line,0] = new_id
        where_to_replace = np.where(prev_angle.T[1] == old_id)
        if len(where_to_replace[0]) > 0:
            for line in where_to_replace[0]:
                new_angle[line,1] = new_id 
        where_to_replace = np.where(prev_angle.T[2] == old_id)
        if len(where_to_replace[0]) > 0:
            for line in where_to_replace[0]:
                new_angle[line,2] = new_id 
silica_info.angles = new_angle

## write silica force field 
f = open('ff/silica.itp','w')

f.write('[ atomtypes]\n')
f.write(';name at.num	mass	charge	ptype	sigma	epsilon\n');
f.write('OC23	8	15.99940	0.000000	A	0.309141855195	0.225936 \n');
f.write('SC4	14	28.08600	0.000000	A	0.369722968028	0.389112\n ')
f.write('HOY	1	1.008000	0.000000	A	0.096662510918	0.062760\n ')
f.write('OC24	8	15.99940	0.000000	A	0.309141855195	0.510448\n\n ')

f.write('[ bondtypes ]\n')
f.write('; i	j	func	b0	kb\n')
f.write('SC4	OC23	1	0.1680	238488.00\n')
f.write('SC4	OC24	1	0.1680	238488.00\n')
f.write('OC24	HOY	1	0.0945	414216.00\n\n')

f.write('[ angletypes ]\n')
f.write('; i	j	k	func	th0	cth\n')
f.write('SC4	OC23	SC4	5	149.00	836.8000	0.0	0.0\n')
f.write('SC4	OC24	HOY	5	115.00	418.4000	0.0	0.0\n')
f.write('OC23	SC4	OC24	5	109.50	836.8000	0.0	0.0\n')
f.write('OC23	SC4	OC23	5	109.50	836.8000	0.0	0.0\n\n')

f.write('[ moleculetype ]\n')
f.write('; molname    nrexcl\n')
f.write('SiOH    3\n')

f.write('\n')
f.write('[ atoms ]\n')
f.write(';   nr       type  resnr residue  atom   cgnr     charge       mass\n')
for cpt in range(len(atom_mass)):
    f.write("{: >6}".format(str(atom_id[cpt])))
    f.write("{: >11}".format(str(atom_type[cpt])))
    f.write("{: >7}".format(str(residue_number[cpt])))
    f.write("{: >7}".format(str(residue_name[cpt])))
    f.write("{: >7}".format(str(atom_name[cpt])))
    f.write("{: >7}".format(str(atom_id[cpt])))
    f.write("{: >11}".format(str(atom_charge[cpt])))
    f.write("{: >11}".format(str(atom_mass[cpt])))
    f.write('\n')

f.write('\n')
f.write('[ bonds ]\n');
f.write(';  ai    aj funct            c0            c1            c2            c3\n');
for cpt, bond in enumerate(silica_info.bonds):
    f.write("{: >5}".format(str(bond[0])))
    f.write("{: >5}".format(str(bond[1])))
    f.write("{: >6}".format(str(1)))
    f.write('\n')

f.write('\n')
f.write('[ angles ]\n');
f.write(';  ai    aj funct            c0            c1            c2            c3\n');
for cpt, angle in enumerate(silica_info.angles):
    f.write("{: >5}".format(str(angle[0])))
    f.write("{: >5}".format(str(angle[1])))
    f.write("{: >5}".format(str(angle[2])))
    f.write("{: >6}".format(str(5)))
    f.write('\n')     

f.write('\n')
f.write('; Include Position restraint file\n')
f.write('#ifdef POSRES\n')
f.write('#include "posre_SiOH.itp"\n')
f.write('#endif\n')
f.close()

## write posres
f = open('ff/posre_SiOH.itp', 'w')
f.write('; position restraints for silica in water\n\n')
f.write('[ position_restraints ]\n')
f.write(';  i funct       fcx        fcy        fcz\n')
for cpt in range(len(atom_mass)):
    f.write("{: >4}".format(str(atom_id[cpt])))
    f.write("{: >5}".format(str(1)))
    f.write("{: >11}".format(str(0)))
    f.write("{: >11}".format(str(0)))
    f.write("{: >11}".format(str(1000)))
    f.write('\n')
f.close()