#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import warnings
warnings.filterwarnings('ignore')


# ### number of segments for the PEG molecules

# In[2]:


Nseg = 8


# ### box size

# In[3]:


Lx, Ly, Lz = 50, 20, 20


# ### load end-groups (head) and monomere up and down

# In[4]:


endpatch = np.loadtxt('DATA/endpatch.dat')
monomer = np.loadtxt('DATA/monomer.dat')
v = 0.28 # distance between 2 monomers


# ### place number of segments

# In[5]:


atoms = np.zeros((10000,7))
cptatoms = 0
# place patch 1
for m in endpatch:
    atoms[cptatoms] = cptatoms+1, 1, m[1], m[2], -m[3]-v, -m[4], m[5]
    cptatoms += 1
# place N monomers
for seg in range(Nseg):
    for m in monomer:
        atoms[cptatoms] = cptatoms+1, 1, m[1], m[2], m[3]+seg*v*2, m[4], m[5]
        cptatoms += 1   
# place patch 2
for m in endpatch:
    atoms[cptatoms] = cptatoms+1, 1, m[1], m[2], m[3]+(2*seg)*v+v, m[4], m[5]
    cptatoms += 1
atoms = atoms[:cptatoms]
atoms.T[4] *= 10 # to Angstroms
atoms.T[5] *= 10 # to Angstroms
atoms.T[6] *= 10 # to Angstroms
car = atoms[atoms.T[2] == 1] # CC32A (bulk carbon)
hyd = atoms[(atoms.T[2] == 3) | (atoms.T[2] == 5)] # 3 = HCA2 (bulk hydrogen), 5 = HCP1 (surface hydrogen)
oxy = atoms[(atoms.T[2] == 2) | (atoms.T[2] == 4)] # 2 = OC30A (bulk hydrogen), 4 = OC311 (surface oxygen)
print(str(len(atoms[(atoms.T[2] == 1)])) + " carbon") 
print(str(len(atoms[(atoms.T[2] == 2)])) + " bulk oxygen")
print(str(len(atoms[(atoms.T[2] == 3)])) + " bulk hydrogen")
print(str(len(atoms[(atoms.T[2] == 4)])) + " surface oxygen")
print(str(len(atoms[(atoms.T[2] == 5)])) + " surface hydrogen")     


# In[6]:


conversion_table = [["CC32A", "OC30A", "HCA2", "OC311", "HCP1"], [1, 2, 3, 4, 5]]


# ### center PEG in box

# In[7]:


atoms.T[4] -= np.mean(atoms.T[4])
atoms.T[5] -= np.mean(atoms.T[5])
atoms.T[6] -= np.mean(atoms.T[6])


# ### estimate molar mass

# In[8]:


molmass = len(car)*12+len(oxy)*16+len(hyd)*1
print('PEG - '+str(molmass)+' g/mol')


# ### add bonds

# In[9]:


bonds = np.zeros((10000,2))
bond_types = np.zeros((10000))
cptbonds = 0
# carbon - carbon bonds between monomers (CC32A - CC32A)
ccbonds = 0
for idx0 in np.int32(car.T[0][:-1:2]):
    idx1 = np.int32(car.T[0][np.where(car.T[0] == idx0)[0][0]+1])
    if idx0<idx1:
        bonds[cptbonds] = idx0, idx1
    else:
        bonds[cptbonds] = idx1, idx0
    bond_types[cptbonds] = 1
    cptbonds += 1
    ccbonds += 1
print(str(ccbonds) + " carbon - carbon bonds")


# In[10]:


# carbon - oxygen bonds
xyz = car.T[4:].T
cobonds = 0
co_surf_bonds = 0
for n0 in range(len(oxy)):
    xyz0 = oxy[n0][4:]
    idx0 = np.int32(oxy[n0][0])
    d = np.sqrt((xyz.T[0]-xyz0[0])**2+(xyz.T[1]-xyz0[1])**2+(xyz.T[2]-xyz0[2])**2)
    where = np.where((d > 0) & (d < 1.5))
    for w in where[0]:
        idx1 = np.int32(car[w][0])
        if idx0<idx1:
            bonds[cptbonds] = idx0, idx1
        else:
            bonds[cptbonds] = idx1, idx0
        type1 =  oxy[n0][2]
        type2 =  car[w][2]
        if (type1 == 4) & (type2 == 1): # surface oxygen with carbon (OC311 - CC32A)
            bond_types[cptbonds] = 4
            co_surf_bonds += 1
        elif (type1 == 2) & (type2 == 1): # bulk oxygen with carbon (OC30A - CC32A)
            bond_types[cptbonds] = 3
            cobonds += 1
        cptbonds += 1
print(str(cobonds) + " bulk oxygen - carbon bond")
print(str(co_surf_bonds) + " surface oxygen - carbon bond")


# In[11]:


# carbon - hydrogen bonds
xyz = car.T[4:].T
chbonds = 0
for n0 in range(len(hyd)):
    xyz0 = hyd[n0][4:]
    idx0 = np.int32(hyd[n0][0])
    d = np.sqrt((xyz.T[0]-xyz0[0])**2+(xyz.T[1]-xyz0[1])**2+(xyz.T[2]-xyz0[2])**2)
    where = np.where((d > 0) & (d < 1.1))[0]
    if where.shape == (1,):
        idx1 = car[where][0][0]
        if idx0<idx1:
            bonds[cptbonds] = idx0, idx1
        else:
            bonds[cptbonds] = idx1, idx0
        type1 =  hyd[n0][2]
        type2 =  car[where[0]][2] 
        assert type1 == 3
        assert type2 == 1
        bond_types[cptbonds] = 3 # bulk carbon with hydrogen (CC32A - HCA2)
        cptbonds += 1
        chbonds += 1
print(str(chbonds) + " bulk hydrogen - carbon bond")


# In[12]:


# oxygen - hydrogen bonds
xyz = oxy.T[4:].T
ohbonds = 0
for n0 in range(len(hyd)):
    xyz0 = hyd[n0][4:]
    idx0 = np.int32(hyd[n0][0])
    d = np.sqrt((xyz.T[0]-xyz0[0])**2+(xyz.T[1]-xyz0[1])**2+(xyz.T[2]-xyz0[2])**2)
    where = np.where((d > 0) & (d < 1.1))[0]
    if where.shape == (1,):
        idx1 = oxy[where][0][0]
        if idx0<idx1:
            bonds[cptbonds] = idx0, idx1
        else:
            bonds[cptbonds] = idx1, idx0
        type1 =  hyd[n0][2]
        type2 =  oxy[where[0]][2] 
        bond_types[cptbonds] = 5 # surface oxygen with hydrogen (OC311 - HCP1)            
        cptbonds += 1  
        ohbonds += 1
#bonds = bonds[bonds[:, 0].argsort()]
print(str(ohbonds) + " surface oxygen - hydrogen bond")


# In[13]:


# remove excess lines and reorder
bonds = bonds[:cptbonds]
bond_types = bond_types[:cptbonds]


# ### calculate angles

# In[14]:


angles = np.zeros((10000,3))
angle_types = np.zeros(10000)
cptangles = 0
bonded_a = np.append(bonds.T[0],bonds.T[1])
for a in atoms:
    ida = np.int32(a[0])
    tpa = np.int32(atoms[atoms.T[0] == ida].T[2])[0]
    occurence = np.sum(bonded_a == ida)
    if occurence > 1: # the atom has 2 or more neighbors
        id_neighbors = np.unique(bonds[(bonds.T[0] == ida) | (bonds.T[1] == ida)].T[:2].T)
        for idb in id_neighbors:
            for idc in id_neighbors:
                if (idb != ida) & (idc != ida) & (idb < idc): # avoid counting same angle twice
                    angles[cptangles] = idb, ida, idc
                    tpb = np.int32(atoms[atoms.T[0] == idb].T[2])[0]
                    tpc = np.int32(atoms[atoms.T[0] == idc].T[2])[0]
                    if (tpb == 3) & (tpa == 1) & (tpc == 3): # HCA2 CC32A HCA2
                        angle_types[cptangles] = 5
                    elif ((tpb == 3) & (tpa == 1) & (tpc == 4)) | ((tpb == 4) & (tpa == 1) & (tpc == 3)): # HCA2 CC32A OC311
                        angle_types[cptangles] = 6
                    elif ((tpb == 3) & (tpa == 1) & (tpc == 1)) | ((tpb == 1) & (tpa == 1) & (tpc == 3)): # HCA2 CC32A CC32A
                        angle_types[cptangles] = 2
                    elif ((tpb == 4) & (tpa == 1) & (tpc == 1)) | ((tpb == 1) & (tpa == 1) & (tpc == 4)): # OC311 CC32A CC32A
                         angle_types[cptangles] = 7
                    elif ((tpb == 1) & (tpa == 4) & (tpc == 5)) | ((tpb == 5) & (tpa == 4) & (tpc == 1)): # CC32A OC311 HCP1    
                         angle_types[cptangles] = 1
                    elif ((tpb == 1) & (tpa == 1) & (tpc == 2)) | ((tpb == 2) & (tpa == 1) & (tpc == 1)): # CC32A CC32A OC30A     
                         angle_types[cptangles] = 3
                    elif ((tpb == 3) & (tpa == 1) & (tpc == 2)) | ((tpb == 2) & (tpa == 1) & (tpc == 3)): # HCA2 CC32A OC30A    
                         angle_types[cptangles] = 4
                    elif (tpb == 1) & (tpa == 2) & (tpc == 1): # CC32A OC30A CC32A
                         angle_types[cptangles] = 8
                    else:   
                        print(tpb, tpa, tpc)
                        print(conversion_table[0][tpb-1], conversion_table[0][tpa-1], conversion_table[0][tpc-1])
                        print("Unknown angle")
                    cptangles += 1       
angles = angles[:cptangles]
angle_types = angle_types[:cptangles]


# ## calculate dihedrals

# In[15]:


dihedrals = np.zeros((10000,4))
dihedral_types = np.zeros(10000)
cptdihedrals = 0
central_angled_a = angles.T[1]
edge_angled_a = np.append(angles.T[0],angles.T[2])
for a in atoms:
    ida = np.int32(a[0])
    tpa = np.int32(atoms[atoms.T[0] == ida].T[2])[0]
    if (tpa == 1) | (tpa == 2) | (tpa == 4): # ignore hydrogen
        id_first_neighbor = np.unique(angles[(angles.T[1] == ida)].T[:3].T)
        id_first_neighbor = id_first_neighbor[id_first_neighbor != ida]
        for idb in id_first_neighbor:
            id_second_neighbor = np.unique(angles[(angles.T[1] == idb)].T[:3].T)
            if len(id_second_neighbor)>0:
                id_second_neighbor = id_second_neighbor[id_second_neighbor != idb]
                id_second_neighbor = id_second_neighbor[id_second_neighbor != ida]
                for idc in id_first_neighbor:
                    if idc != idb:
                        for ide in id_second_neighbor:
                            tpc = np.int32(atoms[atoms.T[0] == idc].T[2])[0]
                            tpb = np.int32(atoms[atoms.T[0] == idb].T[2])[0]
                            tpe = np.int32(atoms[atoms.T[0] == ide].T[2])[0]
                            if (ida < idb) & (tpc != 3) & (tpe != 3) : 
                                if ((tpc == 1) & (tpa == 1) & (tpb == 4) & (tpe == 5)) |                                    ((tpc == 5) & (tpa == 4) & (tpb == 1) & (tpe == 1)):
                                    dihedral_types[cptdihedrals] = 6
                                    dihedrals[cptdihedrals] = idc, ida, idb, ide
                                    cptdihedrals += 1
                                    dihedral_types[cptdihedrals] = 7
                                    dihedrals[cptdihedrals] = idc, ida, idb, ide
                                    cptdihedrals += 1
                                    dihedral_types[cptdihedrals] = 8
                                    dihedrals[cptdihedrals] = idc, ida, idb, ide
                                    cptdihedrals += 1                                    
                                elif ((tpc == 4) & (tpa == 1) & (tpb == 1) & (tpe == 2)) |                                      ((tpc == 2) & (tpa == 1) & (tpb == 1) & (tpe == 4)):
                                    dihedral_types[cptdihedrals] = 12
                                    dihedrals[cptdihedrals] = idc, ida, idb, ide
                                    cptdihedrals += 1
                                    dihedral_types[cptdihedrals] = 13
                                    dihedrals[cptdihedrals] = idc, ida, idb, ide
                                    cptdihedrals += 1
                                    dihedral_types[cptdihedrals] = 14
                                    dihedrals[cptdihedrals] = idc, ida, idb, ide
                                    cptdihedrals += 1                                    
                                elif ((tpc == 1) & (tpa == 1) & (tpb == 2) & (tpe == 1)) |                                      ((tpc == 1) & (tpa == 2) & (tpb == 1) & (tpe == 1)):
                                    dihedral_types[cptdihedrals] = 1
                                    dihedrals[cptdihedrals] = idc, ida, idb, ide
                                    cptdihedrals += 1
                                    dihedral_types[cptdihedrals] = 2
                                    dihedrals[cptdihedrals] = idc, ida, idb, ide
                                    cptdihedrals += 1
                                    dihedral_types[cptdihedrals] = 3
                                    dihedrals[cptdihedrals] = idc, ida, idb, ide
                                    cptdihedrals += 1                                    
                                elif (tpc == 2) & (tpa == 1) & (tpb == 1) & (tpe == 2) :
                                    dihedral_types[cptdihedrals] = 4
                                    dihedrals[cptdihedrals] = idc, ida, idb, ide
                                    cptdihedrals += 1
                                    dihedral_types[cptdihedrals] = 5
                                    dihedrals[cptdihedrals] = idc, ida, idb, ide
                                    cptdihedrals += 1
                                else:                                     
                                    print(tpc, tpa, tpb, tpe)
                                    print(conversion_table[0][tpc-1], conversion_table[0][tpa-1], 
                                          conversion_table[0][tpb-1], conversion_table[0][tpe-1])
                                    print("Unknown angle")     
dihedrals = dihedrals[:cptdihedrals]
dihedral_types = dihedral_types[:cptdihedrals]


# In[16]:


# space for water
bond_types += 1
angle_types += 1
atoms.T[2] += 2


# ## write data.lammps file

# In[17]:


f = open("data.lammps", "w")
f.write('# LAMMPS data file \n\n')
f.write(str(len(atoms))+' atoms\n')
f.write(str(len(bonds))+' bonds\n')
f.write(str(len(angles))+' angles\n')
f.write(str(len(dihedrals))+' dihedrals\n')
f.write('\n')
f.write('7 atom types\n')
f.write('6 bond types\n')
f.write('9 angle types\n')
f.write('14 dihedral types\n')
f.write('\n')
f.write(str(-Lx/2)+' '+str(Lx/2)+' xlo xhi\n')
f.write(str(-Ly/2)+' '+str(Ly/2)+' ylo yhi\n')
f.write(str(-Lz/2)+' '+str(Lz/2)+' zlo zhi\n')
f.write('\n')
f.write('Atoms\n')
f.write('\n')
for myatom in atoms:
    for col in range(len(myatom)):
        if col < 3:
            f.write(str(int(myatom[col]))+' ')
        else :
            f.write(str(myatom[col])+' ')
    f.write('\n')
f.write('\n')
f.write('Bonds\n')
f.write('\n')
for cpt, mybond in enumerate(bonds):
    myline = [cpt + 1, bond_types[cpt], mybond[0], mybond[1]]
    for col in range(len(myline)):
        f.write(str(int(myline[col]))+' ')
    f.write('\n')
f.write('\n')
f.write('Angles\n')
f.write('\n')   
for cpt, myangle in enumerate(angles):
    myline = [cpt + 1, angle_types[cpt], myangle[0], myangle[1], myangle[2]]
    for col in range(len(myline)):
        f.write(str(int(myline[col]))+' ')
    f.write('\n')
f.write('\n')
f.write('Dihedrals\n')
f.write('\n')   
for cpt, mydihedral in enumerate(dihedrals):
    myline = [cpt + 1, dihedral_types[cpt], mydihedral[0], mydihedral[1], mydihedral[2], mydihedral[3]]
    for col in range(len(myline)):
        f.write(str(int(myline[col]))+' ')
    f.write('\n')
f.close()

