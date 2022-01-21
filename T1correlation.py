#!/usr/bin/env python
# coding: utf-8

# In[1]:


import gc
import time
import copy
import random
import os.path
import warnings
import numpy as np
from tqdm import tqdm
import MDAnalysis as mda
from pathlib import Path
from numpy.linalg import norm
from scipy.special import sph_harm
from scipy import constants as cst
from matplotlib import pyplot as plt
import MDAnalysis.lib.distances as dst
from scipy.interpolate import interp1d
warnings.filterwarnings("ignore")


# In[2]:


def FourierTransform(data, method="np.rfft"):
    """
    Wrapper function that takes the data in real space with
    columns time, frequency and returns the data in Fourier space
    with columns frequency, signal
    Units in : ps, signal
    Units out : MHz, s*signal
    """
    dt = (data[1, 0] - data[0, 0]) * cst.pico # second
    return np.vstack(
    (np.fft.rfftfreq(len(data), dt) / cst.mega, # f in MHz
    np.fft.rfft(data[:, 1], norm=None) * dt * 2)).T


# In[3]:


def CorrelationFunction(a,b=None):
    """Uses fast fourier transforms to give the correlation function
    of two arrays, or, if only one array is given, the autocorrelation."""    
    if len(a.shape) > 1:
        a2 = np.append(a, np.zeros((2**int(np.ceil((np.log(len(a))/np.log(2))))-len(a),a.shape[1])),axis=0)
    else:
        a2 = np.append(a, np.zeros(2**int(np.ceil((np.log(len(a))/np.log(2))))-len(a)),axis=0)
    data_a = np.append(a2, np.zeros(a2.shape),axis=0)
    fra = np.fft.fft(data_a,axis=0)
    
    if b is None:
        sf = np.conj(fra)*fra
    else:
        if len(b.shape) > 1:
            b2 = np.append(b, np.zeros((2**int(np.ceil((np.log(len(b))/np.log(2))))-len(b),b.shape[1])),axis=0)
        else:
            b2 = np.append(b, np.zeros(2**int(np.ceil((np.log(len(b))/np.log(2))))-len(b)),axis=0)
        data_b = np.append(b2, np.zeros(b2.shape),axis=0)
        frb = np.fft.fft(data_b,axis=0)
        sf = np.conj(fra)*frb
    res = np.fft.ifft(sf,axis=0)
    
    if len(a.shape) > 1:
        # average over all particles/molecules
        cor = (np.real(res[:len(a)])/np.array(range(len(a),0,-1))[:,np.newaxis]).mean(axis=1)
    else:
        cor = np.real(res[:len(a)])/np.array(range(len(a),0,-1))

    return cor


# In[4]:


def SecondMoment(C):
    """
    Calculate second moment Delta omega
    G(0) = (1/3) (Delta omega)**2
    Unit in : Hz**2
    Unit out : Hz
    """
    if C.shape[0] > 1:
        if len(C.T[0]) == 3:
            Dw = np.sqrt(3*C.T[0])
        elif len(C[0]) == 3:
            Dw = np.sqrt(3*C[0])
        else:
            print('format of the correlation data not understood')
    elif C.shape[0] == 1 :
        Dw = np.sqrt(3*C[0][0])
    else:
        print('format of the correlation data not understood')
    return np.real(Dw)


# In[5]:


def CorrelationTime(C,J,Ct,cut):
    """
    Calculate correlation time
    tau = 0.5 J(0) / G(0)
    tau = 1/G(0) int_0^inf G(t) dt
    Unit in : Hz**2, Hz
    Unit out : s
    """
    if C.shape[0] > 1 :
        if (len(C.T[0]) == 3) & (len(J.T[0]) == 3):
            tau1 = 0.5*(J.T[0]/C.T[0])
            tau2 = np.trapz(C.T[:cut].T,t[:cut]*1e-12)/C.T[0]
        elif (len(C[0]) == 3) & (len(J[0]) == 3):
            tau1 = 0.5*(J[0]/C[0])
            tau2 = np.trapz(C[:cut],t[:cut]*1e-12)/C[0]
        else:
            print('format of C & J data not understood')
    elif C.shape[0] == 1 :
        tau1 = 0.5*(J[0][0]/C[0][0])
        tau2 = np.trapz(C[0][:cut],t[:cut]*1e-12)/C[0][0]
    else:
        print('format of C & J data not understood')
    return np.real(tau1), np.real(tau2)


# In[6]:


def ImportUnivers(fileName="run", mystart=0, myend=0 , mystep=1):
    """import MDA univers"""
    grofile = 'conf.gro'
    trjfile = fileName+'.xtc'
    tprfile = fileName+'.tpr'
    if (os.path.isfile(tprfile)) & (os.path.isfile(trjfile)): # use trp file
        u = mda.Universe(tprfile, trjfile)
    elif (os.path.isfile(grofile)) & (os.path.isfile(trjfile)): # use conf file
        u = mda.Universe(grofile, trjfile)
    else:
        print('trajectory or/and topol file(s) missing')
        return 0
    if (mystart>0) | (myend>0):
        u.transfer_to_memory(start=mystart, stop=myend)
    elif mystep>1:
        u.transfer_to_memory(step=step)
    elif ((mystart>0) | (myend>0)) & (mystep>1):
        u.transfer_to_memory(start=mystart, stop=myend, step=mystep)
    return u


# In[7]:


def SelectProton(u, typeI, group=0):
    """
    select protons from the univers
    if typeI == 'R', grp1 and grp2 are proton from the same water molecules
    if typeI = 'T', grp1 is one single proton, and grp2 is all other atoms not
    in the same molecule
    """
    
    if typeI == 'R':
        resids = ' '.join(str(e) for e in group)
        grp1 = u.select_atoms('name HW1 and resid '+resids)
        grp2 = u.select_atoms('name HW2 and resid '+resids)        
        assert grp1.atoms.n_atoms == grp2.atoms.n_atoms
        assert np.unique(grp1.atoms.names)[0] == 'HW1'
        assert np.unique(grp2.atoms.names)[0] == 'HW2'
        assert np.unique(grp1.atoms.resindices == grp2.atoms.resindices)[0] == True

    elif typeI == 'T':
        grpH = u.select_atoms('type HW')
        idx1 = np.array(random.choices(grpH.atoms.indices,k = 1))
        grp1 = u.select_atoms('index '+str(idx1[0]))
        rH1 = grpH.resids[grpH.atoms.indices == idx1[0]]
        assert grp1.atoms.n_atoms == 1
        idx2 = grpH.atoms.indices[grpH.resids != rH1]
        str2 = ' '.join(str(e) for e in idx2)
        grp2 = u.select_atoms('index '+str2)
        assert grp2.atoms.n_atoms == grpH.atoms.n_atoms-2   
    else:
        print('unknown selection')
        
    assert grp1.atoms.n_atoms>0, 'no atom in group 1'
    assert grp2.atoms.n_atoms>0, 'no atom in group 2'
        
    return grp1, grp2


# In[8]:


def CartesianToSpherical(x,y,z):
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan2(np.sqrt(x**2 + y**2), z)
    phi = np.arctan2(y, x)
    return r, theta, phi


# In[9]:


def SphericalHarmonic(rthephi,order):
    """evaluate harmonic functions
    convention : theta = polar angle, phi = azimuthal angle
    note: scipy uses the opposite convention
    """
    if (order[0]>0) & (order[1] == 0) & (order[2] == 0):
        sh20 = sph_harm(0, 2, rthephi[2], rthephi[1])/np.power(rthephi[0],3)  
        return sh20
    elif np.sum(order) == 3:
        sh20 = sph_harm(0, 2, rthephi[2], rthephi[1])/np.power(rthephi[0],3) 
        sh21 = sph_harm(1, 2, rthephi[2], rthephi[1])/np.power(rthephi[0],3)
        sh22 = sph_harm(2, 2, rthephi[2], rthephi[1])/np.power(rthephi[0],3)
        return sh20, sh21, sh22
    else:
        print('Combination of order not implemented')


# In[10]:


gamma = 2*np.pi*42.6e6 # gyromagnetic constant in Hz/T
K = (3*np.pi/5)*(cst.mu_0/4/np.pi)**2*(cst.hbar)**2*(gamma)**4 # m6/s2
def CorrelationProton(u, grp1, grp2, order,typeI):
    """
    loop on trajectories,
    evaluate distances and angles between atoms,
    calculate harmonic functions,
    then measure correlation functions of harmonic functions
    """
    
    natom = grp2.atoms.n_atoms
    nfrm = u.trajectory.n_frames
    M = np.zeros((np.sum(order),nfrm,natom), dtype=np.complex64)
    C = np.zeros((np.sum(order),nfrm),dtype=np.float32)

    cpt = 0
    for ts in u.trajectory:
        posH1 = grp1.atoms.positions
        posH2 = grp2.atoms.positions
        box = u.dimensions
        # vector between grp1 and grp2 with pbc
        xyz = (np.remainder(posH2 - posH1 + box[:3]/2., box[:3]) - box[:3]/2.).T
        # convert cartesian to sphere
        rthephi = CartesianToSpherical(xyz[0], xyz[1], xyz[2])
        # evaluate spherical harmonic / r**3
        M[:,cpt] = SphericalHarmonic(rthephi,order)
        cpt += 1
        
        if typeI == 'R':
            assert np.max(rthephi[0])<1.8, 'distance too large within molecules'
        assert len(rthephi[0]) == natom, 'weird numer of atoms in group'

    for m in range(np.sum(order)):
        for idx in range(natom):
            C[m] += CorrelationFunction(M[m,:,idx])
        if typeI == 'R':
            C[m] /= natom # normalise by the number of atoms in the group
    C *= K / cst.angstrom**6 # convert to m6/s2 units
    return C


# ### Import univers

# In[16]:


tini = 30000
tfin = 40000
u = ImportUnivers(mystart=tini, myend=tfin)
order_m = [1,1,1] # calculate higher order m=1 and m=2
cut = np.min([u.trajectory.n_frames,1024]) # cut off for correlation function (impact correlation time)
t = np.arange(cut)*np.round(u.trajectory.dt,3) # time vector for correlation function in ps


# ### Intra molecular contribution $R$

# In[29]:


NR = 50 # number of slice for CR calculation (for better performances, does not impact the result)
resOw = u.select_atoms("name OW").atoms.resids # isolate H2O resids ids
groups = np.split(resOw,NR) # split the water molecules in NR groups
CR = np.zeros((np.sum(order_m),u.trajectory.n_frames), dtype=np.float32) # allocate memory
JR = np.zeros((np.sum(order_m),cut//2+1), dtype=np.float32) # allocate memory
for grp in tqdm(groups): # loop on the sub-group
    grp1, grp2 = SelectProton(u, 'R', group=grp)
    CR += CorrelationProton(u, grp1, grp2, order_m,'R')
    gc.collect()
CR /= len(groups)
for m in range(np.sum(order_m)): # evaluate spectrum
    f, JR[m] = np.real(FourierTransform(np.vstack([t[:cut],CR[m][:cut]]).T).T)
DwR = SecondMoment(CR)
tauRa, tauRb = CorrelationTime(CR,JR,t,cut)
print('The second moment (R) is '+str(np.round(DwR/2/np.pi/1000,2))+' kHz')
print('The correla. time (R) is '+str(np.round(tauRa*1e12,3))+' ps')


# ### Inter molecular contribution $T$

# In[30]:


NT = 50 # number of atoms considered for CT calculation
CT = np.zeros((np.sum(order_m),u.trajectory.n_frames), dtype=np.float32) # allocate memory
JT = np.zeros((np.sum(order_m),cut//2+1), dtype=np.float32) # allocate memory
for ids in tqdm(range(NT)):
    grp1, grp2 = SelectProton(u, typeI='T')
    CT += CorrelationProton(u, grp1, grp2, order_m,'T')
    gc.collect()
CT /= NT
for m in range(np.sum(order_m)): # evaluate spectrum
    f, JT[m] = np.real(FourierTransform(np.vstack([t[:cut],CT[m][:cut]]).T).T)
DwT = SecondMoment(CT)
tauTa, tauTb = CorrelationTime(CT,JT,t,cut)
print('The second moment (T) is '+str(np.round(DwT/2/np.pi/1000,2))+' kHz')
print('The correla. time (T) is '+str(np.round(tauTa*1e12,3))+' ps')


# ### Estimate $T_1$

# In[31]:


T1R = 1/(10/3*(DwR)**2*(tauRa))
T1T = 1/(10/3*(DwT)**2*(tauTa))
T1 = 1/(1/T1R+1/T1T)
print('T1 is '+str(np.round(T1,2))+' s')


# ### Calculate $T_1$ and $T_2$ from the spectrums $R_1$ and $R_2$ 
# 
# $ \dfrac{1}{T_1} = J(\omega) + 4 J(2 \omega) $
# 
# $ \dfrac{1}{T_2} = \dfrac{3}{2} J(0) + \dfrac{5}{2} J(\omega) + J(2 \omega) $

# In[33]:


intJR20 = interp1d(f, JR[0], fill_value="extrapolate")
intJT20 = interp1d(f, JT[0], fill_value="extrapolate")
R1R = (1/(intJR20(f) + 4*intJR20(2*f)))
R1T = (1/(intJT20(f) + 4*intJT20(2*f)))
R1 = 1/(1/R1R+1/R1T)

R2R = (1/((3/2)*intJR20(0) + (5/2)*intJR20(f)+intJR20(2*f)))
R2T = (1/((3/2)*intJT20(0) + (5/2)*intJT20(f)+intJT20(2*f)))
R2 = 1/(1/R2R+1/R2T)
#print('T1 is '+str(np.round(R1[0],2))+' s')
#print('T2 is '+str(np.round(R2[0],2))+' s')


# ### Save data into dictionary

# In[36]:


box = u.dimensions
dictionary = {
    "dt": u.trajectory.dt,
    "Lx": box[0],
    "Ly": box[1],
    "Lz": box[2],
    "t": t,
    "f" :f,
    "CR": CR,
    "CT": CT,
    "JR": JR,
    "JT": JT,
    "DwR": DwR,
    "DwT": DwT,
    "tauRa": tauRa,
    "tauTa": tauTa,
    "tauRb": tauRb,
    "tauTb": tauTb,
    "cut": cut,
    "NT": NT,
    "R1":R1,
    "R2":R2,
    "nmol": u.atoms.n_atoms,
    "T1R":T1R, 
    "T1T":T1R, 
    "T1":T1}
np.save('results_correlation_'+str(tini)+'_'+str(tfin)+'.npy', dictionary)


# In[ ]:


# Uncomment to save multiple files with different names
#chosen = False
#Ns = 0
#while chosen == False:
#path = 'results'+str(Ns)+'.npy'
#myfile = Path(path)
#if myfile.is_file() == False:
#    np.save(path, dictionary)
#    chosen = True
#else:
#    Ns += 1

