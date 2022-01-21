import numpy as np

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

def CartesianToSpherical(x,y,z):
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan2(np.sqrt(x**2 + y**2), z)
    phi = np.arctan2(y, x)
    return r, theta, phi

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
