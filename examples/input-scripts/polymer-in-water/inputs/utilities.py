import numpy as np
import random
import copy
from numpy.linalg import norm
from scipy.spatial.transform import Rotation as R

def PEGgenerator(Nseg):
    """Generate PEG molecule with desired number of segments."""
    ### load end-groups (head) and monomere up and down
    endpatch = np.loadtxt('PEG/endpatch.dat')
    monomer = np.loadtxt('PEG/monomer.dat')
    v = 0.28 # distance between 2 monomers
    ### place number of segments
    atoms = np.zeros((10000,6))
    cptatoms = 0
    # place patch 1
    for m in endpatch:
        atoms[cptatoms] = cptatoms+1, m[1], m[2], -m[3]-v, -m[4], m[5]
        cptatoms += 1
    # place N monomers
    for seg in range(Nseg):
        for m in monomer:
            atoms[cptatoms] = cptatoms+1, m[1], m[2], m[3]+seg*v*2, m[4], m[5]
            cptatoms += 1   
    # place patch 2
    for m in endpatch:
        atoms[cptatoms] = cptatoms+1, m[1], m[2], m[3]+(2*seg)*v+v, m[4], m[5]
        cptatoms += 1
    atoms = atoms[:cptatoms]
    car = atoms[atoms.T[1] == 1]
    hyd = atoms[(atoms.T[1] == 3) | (atoms.T[1] == 5)]
    oxy = atoms[(atoms.T[1] == 2) | (atoms.T[1] == 4)]
    ### add bonds
    bonds = np.zeros((10000,2))
    cptbonds = 0
    # carbon - carbon bonds between monomers
    for idx0 in np.int32(car.T[0][:-1:2]):
        idx1 = np.int32(car.T[0][np.where(car.T[0] == idx0)[0][0]+1])
        if idx0<idx1:
            bonds[cptbonds] = idx0, idx1
        else:
            bonds[cptbonds] = idx1, idx0
        cptbonds += 1
    # carbon - oxygen bonds
    xyz = car.T[3:].T
    for n0 in range(len(oxy)):
        xyz0 = oxy[n0][3:]
        idx0 = np.int32(oxy[n0][0])
        d = np.sqrt((xyz.T[0]-xyz0[0])**2+(xyz.T[1]-xyz0[1])**2+(xyz.T[2]-xyz0[2])**2)
        where = np.where((d > 0) & (d < 0.15))
        for w in where[0]:
            idx1 = np.int32(car[w][0])
            if idx0<idx1:
                bonds[cptbonds] = idx0, idx1
            else:
                bonds[cptbonds] = idx1, idx0
            cptbonds += 1
    # carbon - hydrogen bonds
    xyz = car.T[3:].T
    for n0 in range(len(hyd)):
        xyz0 = hyd[n0][3:]
        idx0 = np.int32(hyd[n0][0])
        d = np.sqrt((xyz.T[0]-xyz0[0])**2+(xyz.T[1]-xyz0[1])**2+(xyz.T[2]-xyz0[2])**2)
        where = np.where((d > 0) & (d < 0.11))[0]
        if where.shape == (1,):
            idx1 = car[where][0][0]
            if idx0<idx1:
                bonds[cptbonds] = idx0, idx1
            else:
                bonds[cptbonds] = idx1, idx0
            cptbonds += 1      
    # oxygen - hydrogen bonds
    xyz = oxy.T[3:].T
    for n0 in range(len(hyd)):
        xyz0 = hyd[n0][3:]
        idx0 = np.int32(hyd[n0][0])
        d = np.sqrt((xyz.T[0]-xyz0[0])**2+(xyz.T[1]-xyz0[1])**2+(xyz.T[2]-xyz0[2])**2)
        where = np.where((d > 0) & (d < 0.11))[0]
        if where.shape == (1,):
            idx1 = oxy[where][0][0]
            if idx0<idx1:
                bonds[cptbonds] = idx0, idx1
            else:
                bonds[cptbonds] = idx1, idx0
            cptbonds += 1       
    # remove excess lines and reorder
    bonds = bonds[:cptbonds]
    bonds = bonds[bonds[:, 0].argsort()]
    ### calculate angles
    angles = np.zeros((10000,3))
    cptangles = 0
    bonded_a = np.append(bonds.T[0],bonds.T[1])
    for a in atoms:
        ida = np.int32(a[0])
        tpa = np.int32(atoms[atoms.T[0] == ida].T[1])[0]
        occurence = np.sum(bonded_a == ida)
        if occurence > 1: # the atom has 2 or more atoms
            id_neighbors = np.unique(bonds[(bonds.T[0] == ida) | (bonds.T[1] == ida)].T[:2].T)
            for idb in id_neighbors:
                for idc in id_neighbors:
                    if (idb != ida) & (idc != ida) & (idb < idc): # avoid counting same angle twice
                        angles[cptangles] = idb, ida, idc
                        cptangles += 1       
    angles = angles[:cptangles]
    ## calculate dihedrals
    dihedrals = np.zeros((10000,4))
    cptdihedrals = 0
    for a in atoms:
        ida = np.int32(a[0])
        tpa = np.int32(atoms[atoms.T[0] == ida].T[1])[0]
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
                                tpc = np.int32(atoms[atoms.T[0] == idc].T[1])[0]
                                tpe = np.int32(atoms[atoms.T[0] == ide].T[1])[0]
                                if (ida < idb) & (tpc != 3) & (tpe != 3) : 
                                    dihedrals[cptdihedrals] = idc, ida, idb, ide
                                    cptdihedrals += 1
    dihedrals = dihedrals[:cptdihedrals]
    return atoms, bonds, angles, dihedrals

def neighborsearch(neighbor, molecule, x, y, z, Lx, Ly, Lz):
    """Search all neighbor to a molecule in a box and return the closest distance."""
    box = np.array([Lx, Ly, Lz])
    minr = 10
    for m in molecule:
        x0 = m[0] + x
        y0 = m[1] + y
        z0 = m[2] + z
        dxdydz = np.remainder(neighbor - np.array([x0,y0,z0]) + box/2., box) - box/2.
        minr = np.min([minr,np.min(norm(dxdydz,axis=1))])
    return minr

def randomlocation(Lx,Ly,Lz):
    """Choose a random location within a given box."""
    txlo, txhi = -Lx/2, Lx/2
    tylo, tyhi = -Ly/2, Ly/2
    tzlo, tzhi = -Lz/2, Lz/2    
    x = random.randint(1,1000)/1000*(txhi-txlo)
    y = random.randint(1,1000)/1000*(tyhi-tylo)
    z = random.randint(1,1000)/1000*(tzhi-tzlo)
    return x, y, z

def place_molecules(Np, NH2O, atomsPEG, bondsPEG, anglesPEG, dihedralsPEG):
    attempt = 0
    cptPEG = 0
    cptH2O = 0
    while (cptPEG<Np) | (cptH2O<NH2O):
        ## choose initial box dimensions
        Lx, Ly, Lz = 2.5+0.5*attempt, 2.5+0.5*attempt, 2.5+0.5*attempt # nm
        
        box = np.array([Lx, Ly, Lz])
        ## initialise matrix
        atoms = np.zeros((100000,7))
        bonds = np.zeros((100000,2))
        angles = np.zeros((100000,3))
        dihedrals = np.zeros((100000,4))
        cptatoms = 0
        cptbonds = 0
        cptangles = 0
        cptdihedrals = 0
        cptres = 0
        cptPEG = 0
        resName = ["" for x in range(1000000)]
        atoName = ["" for x in range(1000000)]

        fail_attempt = 0
        while (cptPEG < Np) & (fail_attempt<1e3):

            nO = 0
            nH = 0
            nC = 0

            insert = 1
            x,y,z = randomlocation(Lx,Ly,Lz)
            if cptPEG > 0:
                d = neighborsearch(atoms[:cptatoms].T[4:].T,atomsPEG.T[3:].T, x, y, z, Lx, Ly, Lz)
                if d < 1:
                    insert = 0

            if insert == 1:

                for m in bondsPEG:
                    bonds[cptbonds] = m[0]+cptatoms, m[1]+cptatoms
                    cptbonds += 1

                for m in anglesPEG:
                    angles[cptangles] = m[0]+cptatoms, m[1]+cptatoms, m[2]+cptatoms
                    cptangles += 1

                for m in dihedralsPEG:
                    dihedrals[cptdihedrals] = m[0]+cptatoms, m[1]+cptatoms, m[2]+cptatoms, m[3]+cptatoms
                    cptdihedrals += 1

                for m in atomsPEG:
                    atoms[cptatoms] = cptatoms+1, cptres+1, m[1], m[2], m[3]+x, m[4]+y, m[5]+z 
                    resName[cptatoms] = 'PEG'
                    if m[1] == 1:
                        nC += 1
                        atoName[cptatoms] = 'C'+str(nC)
                    elif m[1] == 2:
                        nO += 1
                        atoName[cptatoms] = 'O'+str(nO)              
                    elif m[1] == 3:
                        nH += 1
                        atoName[cptatoms] = 'H'+str(nH)
                    elif m[1] == 4:
                        nO += 1
                        atoName[cptatoms] = 'O'+str(nO)
                    elif m[1] == 5:
                        nH += 1
                        atoName[cptatoms] = 'H'+str(nH)
                    cptatoms += 1
                cptres += 1
                cptPEG += 1

                if cptPEG == 1: ### write PEG.itp

                    f = open('ff/peg.itp', 'w')
                    f.write('[ moleculetype ]\n')
                    f.write('PEG   2\n\n')
                    f.write('[ atoms ]\n')
                    nc = 0
                    no = 0
                    nh = 0
                    for n in range(cptatoms):
                        f.write("{: >5}".format(str(n+1))) # atom number
                        if atoms.T[2][n] == 1:
                            f.write("{: >8}".format('CC32A'))
                        elif atoms.T[2][n] == 2:
                            f.write("{: >8}".format('OC30A'))
                        elif atoms.T[2][n] == 3:
                            f.write("{: >8}".format('HCA2'))
                        elif atoms.T[2][n] == 4:
                            f.write("{: >8}".format('OC311'))
                        elif atoms.T[2][n] == 5:
                            f.write("{: >8}".format('HCP1'))
                        else:
                            print('extra atoms')    
                        f.write("{: >8}".format(str(np.int32(atoms[n][1])))) # residue number
                        f.write("{: >8}".format('PEG')) # residue name
                        if atoms.T[2][n] == 1:
                            nc += 1
                            f.write("{: >8}".format('C'+str(nc))) # atom name
                        elif (atoms.T[2][n] == 3) | (atoms.T[2][n] == 5):
                            nh += 1
                            f.write("{: >8}".format('H'+str(nh))) # atom name
                        elif (atoms.T[2][n] == 2) | (atoms.T[2][n] == 4):
                            no += 1
                            f.write("{: >8}".format('O'+str(no))) # atom name
                        f.write("{: >8}".format(str(np.int32(n+1))))
                        f.write("{: >8}".format(str("{:.3f}".format(atoms.T[3][n]))))
                        if atoms.T[2][n] == 1:
                            f.write("{: >8}".format(str("{:.3f}".format(12.011))))
                        elif (atoms.T[2][n] == 3) | (atoms.T[2][n] == 5):
                            f.write("{: >8}".format(str("{:.3f}".format(1.008))))    
                        elif (atoms.T[2][n] == 2) | (atoms.T[2][n] == 4):
                            f.write("{: >8}".format(str("{:.3f}".format(15.9994)))) 
                        f.write("\n") 
                    f.write("\n")  
                    f.write('[ bonds ]\n')  
                    for n in range(cptbonds):
                        f.write("{: >5}".format(str(np.int32(bonds[n][0]))))
                        f.write("{: >5}".format(str(np.int32(bonds[n][1]))))
                        f.write("{: >5}".format(str(np.int32(1))))
                        f.write("\n")
                    f.write("\n")  
                    f.write('[ angles ]\n')  
                    for n in range(cptangles):
                        f.write("{: >5}".format(str(np.int32(angles[n][0]))))
                        f.write("{: >5}".format(str(np.int32(angles[n][1]))))
                        f.write("{: >5}".format(str(np.int32(angles[n][2]))))
                        f.write("{: >5}".format(str(np.int32(5))))
                        f.write("\n")
                    f.write("\n")  
                    f.write('[ dihedrals ]\n')  
                    for n in range(cptdihedrals):
                        f.write("{: >5}".format(str(np.int32(dihedrals[n][0]))))
                        f.write("{: >5}".format(str(np.int32(dihedrals[n][1]))))
                        f.write("{: >5}".format(str(np.int32(dihedrals[n][2]))))
                        f.write("{: >5}".format(str(np.int32(dihedrals[n][3]))))
                        f.write("{: >5}".format(str(np.int32(9))))
                        f.write("\n")
                    f.close()
            else:
                fail_attempt += 1

        ### place water molecules

        rOH = 0.9572/10 
        rOM = 0.105/10 # tip4p/epsilon
        thetaHOH = 104.52    
        atomH2O = np.array([[1, 6, 0, 0, 0, 0], \
        [2, 7, 0.527, rOH*np.cos((thetaHOH/2)*np.pi/180),   rOH*np.sin((thetaHOH/2)*np.pi/180), 0.0], \
        [3, 7, 0.527, rOH*np.cos((thetaHOH/2)*np.pi/180),   -rOH*np.sin((thetaHOH/2)*np.pi/180),  0.0], \
        [4, 8, -1.054, rOM,  0.0, 0.0]])
        cptH2O = 0
        cptH = 0
        fail_attempt = 0
        while (cptH2O < NH2O) & (fail_attempt <1e6):
            x,y,z = randomlocation(Lx,Ly,Lz)
            pos = np.array([x,y,z])
            dxdydz = np.remainder((pos - atoms[:cptatoms].T[4:].T) + box/2., box) - box/2.
            d = np.min(norm(dxdydz,axis=1))
            if d > 0.31:
                for m in atomH2O:
                    atoms[cptatoms] = cptatoms+1, cptres+1, m[1], m[2], m[3]+x, m[4]+y, m[5]+z
                    resName[cptatoms] = 'SOL'
                    if m[1] == 6:
                        atoName[cptatoms] = 'OW'
                    elif m[1] == 7:
                        if cptH % 2 == 0:
                            atoName[cptatoms] = 'HW1'
                        else:
                            atoName[cptatoms] = 'HW2'
                        cptH += 1
                    elif m[1] == 8:
                        atoName[cptatoms] = 'MW'
                    cptatoms += 1    
                cptH2O += 1
                cptres += 1
            else:
                fail_attempt += 1
        attempt += 1
    
    
    print('box size = '+str(Lx)+' nm')
    atoms = atoms[:cptatoms]
    bonds = bonds[:cptbonds]
    angles = angles[:cptangles]
    dihedrals = dihedrals[:cptdihedrals]
    return atoms, bonds, angles, dihedrals, atoName, resName, Lx, Ly, Lz

def write_topol(Number_polymer, Number_water):

    f = open('topol.top', 'w')
    f.write('#include "ff/forcefield.itp"\n')
    f.write('#include "ff/peg.itp"\n')
    f.write('#include "ff/tip4peps.itp"\n\n')
    f.write('[ System ]\n')
    f.write('PEG in water\n\n')
    f.write('[ Molecules ]\n\n')
    f.write('PEG '+ str(Number_polymer)+'\n')
    f.write('SOL '+ str(Number_water)+'\n')
    f.close()

def write_conf(atoms, atoName, resName, Lx, Ly, Lz):

    f = open('conf.gro', 'w')
    f.write('PEG SYSTEM\n')
    f.write(str(len(atoms))+'\n')
    nc, no, nh, nm = 0,0,0,0
    for n in range(len(atoms)): 
        f.write("{: >5}".format(np.int32(atoms[n][1]))) # residue number (5 positions, integer) 
        f.write("{: >5}".format(str(resName[n]))) # residue name (5 characters) 
        f.write("{: >5}".format(str(atoName[n]))) # residue name (5 characters) 
        f.write("{: >5}".format(str(np.int32(n+1)))) # atom number (5 positions, integer)
        f.write("{: >8}".format(str("{:.3f}".format(atoms[n][4])))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
        f.write("{: >8}".format(str("{:.3f}".format(atoms[n][5])))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places) 
        f.write("{: >8}".format(str("{:.3f}".format(atoms[n][6])))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places) 
        f.write("\n")
    f.write("{: >10}".format(str("{:.5f}".format(Lx))))
    f.write("{: >10}".format(str("{:.5f}".format(Ly))))
    f.write("{: >10}".format(str("{:.5f}".format(Lz))))
    f.write("\n")
    f.close()

def detect_bond_type(t1, t2):
        bond_types = None
        if (t1 == 1) & (t2 == 1): # carbon-carbon (CC32A - CC32A) - 0.15300000 186188.00
            bond_types = 1
        elif ((t1 == 4) & (t2 == 1)) | (t1 == 1) & (t2 == 4): # surface oxygen with carbon (OC311 - CC32A) - 0.14200000 358150.40   
            bond_types = 4
        elif ((t1 == 2) & (t2 == 1)) | ((t1 == 1) & (t2 == 2)): # bulk oxygen with carbon (OC30A - CC32A) - 0.14150000 301248.00
            bond_types = 3
        elif ((t1 == 1) & (t2 == 3)) | ((t1 == 3) & (t2 == 1)): # bulk carbon with hydrogen (CC32A - HCA2) - 0.11110000 258571.20
            bond_types = 2
        elif ((t1 == 4) & (t2 == 5)) | ((t1 == 5) & (t2 == 4)): # surface oxygen with hydrogen (OC311 - HCP1) - 0.09600000 456056.00
            bond_types = 5
        elif ((t1 == 6) & (t2 == 7)) | ((t1 == 7) & (t2 == 6)): # water molecule
            bond_types = 6
        else:
            print("unfound bond")
        return bond_types

conversion_table = [["CC32A", "OC30A", "HCA2", "OC311", "HCP1"], [1, 2, 3, 4, 5]]

def detect_angle_type(tpb, tpa, tpc):
    angle_types = None
    if (tpb == 3) & (tpa == 1) & (tpc == 3): # HCA2 CC32A HCA2
        angle_types = 5
    elif ((tpb == 3) & (tpa == 1) & (tpc == 4)) | ((tpb == 4) & (tpa == 1) & (tpc == 3)): # HCA2 CC32A OC311
        angle_types = 6
    elif ((tpb == 3) & (tpa == 1) & (tpc == 1)) | ((tpb == 1) & (tpa == 1) & (tpc == 3)): # HCA2 CC32A CC32A
        angle_types = 2
    elif ((tpb == 4) & (tpa == 1) & (tpc == 1)) | ((tpb == 1) & (tpa == 1) & (tpc == 4)): # OC311 CC32A CC32A
            angle_types = 7
    elif ((tpb == 1) & (tpa == 4) & (tpc == 5)) | ((tpb == 5) & (tpa == 4) & (tpc == 1)): # CC32A OC311 HCP1    
            angle_types = 1
    elif ((tpb == 1) & (tpa == 1) & (tpc == 2)) | ((tpb == 2) & (tpa == 1) & (tpc == 1)): # CC32A CC32A OC30A     
            angle_types = 3
    elif ((tpb == 3) & (tpa == 1) & (tpc == 2)) | ((tpb == 2) & (tpa == 1) & (tpc == 3)): # HCA2 CC32A OC30A    
            angle_types = 4
    elif (tpb == 1) & (tpa == 2) & (tpc == 1): # CC32A OC30A CC32A
            angle_types = 8
    elif (tpb == 7) & (tpa == 6) & (tpc == 7): # water
            angle_types = 9
    else:   
        print(tpb, tpa, tpc)
        print(conversion_table[0][tpb-1], conversion_table[0][tpa-1], conversion_table[0][tpc-1])
        print("Unknown angle")
    return angle_types

def detect_dihedral_type(tpc, tpa, tpb, tpe, idc, ida, idb, ide):
    dihedral_type = []
    if (ida < idb) & (tpc != 3) & (tpe != 3) : 
        if ((tpc == 1) & (tpa == 1) & (tpb == 4) & (tpe == 5)) | \
            ((tpc == 5) & (tpa == 4) & (tpb == 1) & (tpe == 1)):
            dihedral_type.append(6)
            dihedral_type.append(7)
            dihedral_type.append(8)     
        elif ((tpc == 4) & (tpa == 1) & (tpb == 1) & (tpe == 2)) | \
                ((tpc == 2) & (tpa == 1) & (tpb == 1) & (tpe == 4)):
            dihedral_type.append(12)
            dihedral_type.append(13)
            dihedral_type.append(14)     
        elif ((tpc == 1) & (tpa == 1) & (tpb == 2) & (tpe == 1)) | \
                ((tpc == 1) & (tpa == 2) & (tpb == 1) & (tpe == 1)):
            dihedral_type.append(1)
            dihedral_type.append(2)
            dihedral_type.append(3)             
        elif (tpc == 2) & (tpa == 1) & (tpb == 1) & (tpe == 2) :
            dihedral_type.append(4)
            dihedral_type.append(5)
        else:                                     
            print(tpc, tpa, tpb, tpe)
            print(conversion_table[0][tpc-1], conversion_table[0][tpa-1], 
                    conversion_table[0][tpb-1], conversion_table[0][tpe-1])
            print("Unknown angle") 
    return dihedral_type    

def prepare_lammps(atoms, bonds, angles, dihedrals):
    # GROMACS to LAMMPS
    # conversion_table = [["CC32A", "OC30A", "HCA2", "OC311", "HCP1"], [1, 2, 3, 4, 5]]

    # suppress the extra oxygen
    atoms = atoms[atoms.T[2] < 8]

    # renumber the matrix
    cpt1 = 1
    new_atoms = []
    new_bonds = copy.deepcopy(bonds)
    new_angles = copy.deepcopy(angles)
    new_dihedrals = copy.deepcopy(dihedrals)
    for atom in atoms:
        oldid = atom[0]
        newid = cpt1
        cpt2 = 0
        for bond in bonds:
            if bond[0] == oldid:
                new_bonds[cpt2, 0] = newid
            if bond[1] == oldid:
                new_bonds[cpt2, 1] = newid
            cpt2 += 1
        cpt2 = 0
        for angle in angles:
            if angle[0] == oldid:
                new_angles[cpt2, 0] = newid
            if angle[1] == oldid:
                new_angles[cpt2, 1] = newid
            if angle[2] == oldid:
                new_angles[cpt2, 2] = newid
            cpt2 += 1
        cpt2 = 0
        for dihedral in dihedrals:
            if dihedral[0] == oldid:
                new_dihedrals[cpt2, 0] = newid
            if dihedral[1] == oldid:
                new_dihedrals[cpt2, 1] = newid
            if dihedral[2] == oldid:
                new_dihedrals[cpt2, 2] = newid
            if dihedral[3] == oldid:
                new_dihedrals[cpt2, 3] = newid
            cpt2 += 1
        atom[0] = newid
        # replace 
        new_atoms.append(atom)
        cpt1 += 1
    atoms = np.array(new_atoms)
    bonds = np.array(new_bonds)
    angles = np.array(new_angles)
    dihedrals = np.array(new_dihedrals)
    # add bond for water
    new_bonds = []
    new_angles = []
    cpt = 1
    for atom in atoms:
        if atom[2] == 6:
            new_bonds.append([cpt, cpt+1])
            new_bonds.append([cpt, cpt+2])
            new_angles.append([cpt+1, cpt, cpt+2])
        cpt += 1
    bonds = np.vstack([bonds, new_bonds])
    angles = np.vstack([angles, new_angles])
    return atoms, bonds, angles, dihedrals

def write_lammps(atoms, bonds, angles, dihedrals, Lx, Ly, Lz):

    # estimate dihedral number 
    n_dihedral = 0
    for cpt, mydihedral in enumerate(dihedrals):
        id1 = np.int32(mydihedral[0])
        id2 = np.int32(mydihedral[1])
        id3 = np.int32(mydihedral[2])
        id4 = np.int32(mydihedral[3])
        t1 = np.int32(atoms[atoms.T[0] == id1].T[2][0])
        t2 = np.int32(atoms[atoms.T[0] == id2].T[2][0])
        t3 = np.int32(atoms[atoms.T[0] == id3].T[2][0])
        t4 = np.int32(atoms[atoms.T[0] == id4].T[2][0])
        dihedral_types = detect_dihedral_type(t1, t2, t3, t4, id1, id2, id3, id4)
        for dihedral_type in dihedral_types:
            n_dihedral += 1

    f = open("data.lammps", "w")
    f.write('# LAMMPS data file \n\n')
    f.write(str(len(atoms))+' atoms\n')
    f.write(str(len(bonds))+' bonds\n')
    f.write(str(len(angles))+' angles\n')
    f.write(str(n_dihedral)+' dihedrals\n')
    f.write('\n')
    f.write('7 atom types\n')
    f.write('6 bond types\n')
    f.write('9 angle types\n')
    f.write('14 dihedral types\n')
    f.write('\n')
    f.write(str(-Lx/2*10)+' '+str(Lx/2*10)+' xlo xhi\n')
    f.write(str(-Ly/2*10)+' '+str(Ly/2*10)+' ylo yhi\n')
    f.write(str(-Lz/2*10)+' '+str(Lz/2*10)+' zlo zhi\n')
    f.write('\n')
    f.write('Atoms\n')
    f.write('\n')
    for myatom in atoms:
        if myatom[2] < 8:
            if myatom[2] == 6:
                myatom[3] = -2*0.527
            myatom[4] *= 10
            myatom[5] *= 10
            myatom[6] *= 10
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
        id1 = mybond[0]
        id2 = mybond[1]
        t1 = np.int32(atoms[atoms.T[0] == id1].T[2][0])
        t2 = np.int32(atoms[atoms.T[0] == id2].T[2][0])
        bond_types = detect_bond_type(t1, t2)
        myline = [cpt + 1, bond_types, id1, id2]
        for col in range(len(myline)):
            f.write(str(int(myline[col]))+' ')
        f.write('\n')
    f.write('\n')
    f.write('Angles\n')
    f.write('\n')   
    for cpt, myangle in enumerate(angles):
        id1 = np.int32(myangle[0])
        id2 = np.int32(myangle[1])
        id3 = np.int32(myangle[2])
        t1 = np.int32(atoms[atoms.T[0] == id1].T[2][0])
        t2 = np.int32(atoms[atoms.T[0] == id2].T[2][0])
        t3 = np.int32(atoms[atoms.T[0] == id3].T[2][0])
        angle_types = detect_angle_type(t1, t2, t3)
        if angle_types is not None:
            myline = [cpt + 1, angle_types, id1, id2, id3]
            for col in range(len(myline)):
                f.write(str(int(myline[col]))+' ')
            f.write('\n')
    f.write('\n')
    f.write('Dihedrals\n')
    f.write('\n')   
    n_dihedral = 0
    for cpt, mydihedral in enumerate(dihedrals):
        id1 = np.int32(mydihedral[0])
        id2 = np.int32(mydihedral[1])
        id3 = np.int32(mydihedral[2])
        id4 = np.int32(mydihedral[3])
        t1 = np.int32(atoms[atoms.T[0] == id1].T[2][0])
        t2 = np.int32(atoms[atoms.T[0] == id2].T[2][0])
        t3 = np.int32(atoms[atoms.T[0] == id3].T[2][0])
        t4 = np.int32(atoms[atoms.T[0] == id4].T[2][0])
        dihedral_types = detect_dihedral_type(t1, t2, t3, t4, id1, id2, id3, id4)
        for dihedral_type in dihedral_types:
            myline = [n_dihedral + 1, dihedral_type, mydihedral[0], mydihedral[1], mydihedral[2], mydihedral[3]]
            n_dihedral += 1
            for col in range(len(myline)):
                f.write(str(int(myline[col]))+' ')
            f.write('\n')
    f.close()