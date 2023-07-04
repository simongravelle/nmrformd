#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import random
import copy

Na = 6.022e23 #constants.Avogadro
Mh2o = 0.018053 # kg/mol - water
nwater = 2180
dw = 3.1

# ### loop
attemps = 0
cptH2O = 0
while cptH2O < nwater:
    
    Lx, Ly, Lz = 10+attemps*dw, 10+attemps*dw, 10+attemps*dw
    txlo, txhi = 0, Lx
    tylo, tyhi = 0, Ly
    tzlo, tzhi = 0, Lz
    
    cptatom = 0
    cptbond = 0
    cptangle = 0
    cptmol = 0
    cptH2O = 0
    cptres = 1

    # allocate memory
    XYZ = np.zeros((1000000,3))
    Typ = ["" for x in range(1000000)]
    ResName = ["" for x in range(1000000)]
    ResNum = np.zeros((1000000,1))

    cptH2O = 0
    # create water    
    rOH = 0.9572 
    rOM = 0.1546
    thetaHOH = 104.52    
    PosH2O = np.array([[0, 0, 0],
                       [rOH*np.cos((thetaHOH/2)*np.pi/180), rOH*np.sin((thetaHOH/2)*np.pi/180), 0.0],
                       [rOH*np.cos((thetaHOH/2)*np.pi/180), -rOH*np.sin((thetaHOH/2)*np.pi/180), 0.0],
                       [rOM,  0.0, 0.0]])

    TypH2O = ['OW', 'HW1', 'HW2', 'MW']
    for x in np.arange(txlo+dw/2,txhi,dw):
        for y in np.arange(tylo+dw/2,tyhi,dw):
            for z in np.arange(tzlo+dw/2,tzhi,dw):
                if cptH2O < nwater:
                    for j in range(4):
                        XYZ[cptatom] = [x,y,z]+np.array(PosH2O[j])
                        Typ[cptatom] = TypH2O[j]
                        ResNum[cptatom] = cptres
                        ResName[cptatom] = 'SOL'
                        cptatom += 1    
                    cptH2O += 1
                    cptres += 1
    attemps += 1

print('Lx = '+str(Lx/10)+' nm, Ly = '+str(Ly/10)+' nm, Lz = '+str(Lz/10)+' nm')
print(str(cptH2O)+' water molecules')


# ### write conf.gro
f = open('conf.gro', 'w')
f.write('Pure water\n')
f.write(str(cptatom)+'\n')
for n in range(cptatom):
    f.write("{: >5}".format(str(np.int32(ResNum[n][0])))) # residue number (5 positions, integer) 
    f.write("{: >5}".format(str(ResName[n]))) # residue name (5 characters) 
    f.write("{: >5}".format(str(Typ[n]))) # atom name (5 characters) 
    f.write("{: >5}".format(str(np.int32(n+1)))) # atom number (5 positions, integer)
    f.write("{: >8}".format(str("{:.3f}".format(XYZ[n][0]/10)))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    f.write("{: >8}".format(str("{:.3f}".format(XYZ[n][1]/10)))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places) 
    f.write("{: >8}".format(str("{:.3f}".format(XYZ[n][2]/10)))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places) 
    f.write("\n")
f.write("{: >10}".format(str("{:.5f}".format(Lx/10))))
f.write("{: >10}".format(str("{:.5f}".format(Ly/10))))
f.write("{: >10}".format(str("{:.5f}".format(Lz/10))))
f.close()


# ### write topol.top
f = open('topol.top', 'w')
f.write('#include "ff/forcefield.itp"\n')
f.write('#include "ff/tip4peps.itp"\n')
f.write('[ System ]\n')
f.write('Pure water\n\n')
f.write('[ Molecules ]\n')
f.write('SOL '+ str(cptH2O)+'\n')
f.close()

