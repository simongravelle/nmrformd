 #!/usr/bin/env python3

from scipy import constants as cst
import numpy as np


def fourier_transform(a):
    """
    Wrapper function that takes the data in real space with
    columns time, frequency and returns the data in Fourier space
    with columns frequency, signal
    Units in : ps, signal
    Units out : MHz, s*signal

    Credit to the fourier transform function of MAICoS
    https://maicos-devel.gitlab.io/maicos/index.html
    """
    dt = (a[1, 0] - a[0, 0]) * cst.pico  # second
    return np.vstack((np.fft.rfftfreq(len(a), dt) / cst.mega, np.fft.rfft(a[:, 1], norm=None) * dt * 2)).T


def correlation_function(a, b=None):
    """Uses fast fourier transforms to calculate the correlation function.

    If two arrays are given, the cross-correlation is returned.
    If one array is given, the auto-correlation is returned.

    Credit to the correlation function of MAICoS
    https://maicos-devel.gitlab.io/maicos/index.html
    """
    if len(a.shape) > 1:
        a2 = np.append(a,
                       np.zeros((2 ** int(np.ceil((np.log(len(a)) / np.log(2)))) - len(a), a.shape[1])),
                       axis=0)
    else:
        a2 = np.append(a, np.zeros(2 ** int(np.ceil((np.log(len(a)) / np.log(2)))) - len(a)), axis=0)
    data_a = np.append(a2, np.zeros(a2.shape), axis=0)
    fra = np.fft.fft(data_a, axis=0)
    if b is None:
        sf = np.conj(fra) * fra
    else:
        if len(b.shape) > 1:
            b2 = np.append(b, np.zeros(
                (2 ** int(np.ceil((np.log(len(b)) / np.log(2)))) - len(b), b.shape[1])), axis=0)
        else:
            b2 = np.append(b, np.zeros(2 ** int(np.ceil((np.log(len(b)) / np.log(2)))) - len(b)),
                           axis=0)
        data_b = np.append(b2, np.zeros(b2.shape), axis=0)
        frb = np.fft.fft(data_b, axis=0)
        sf = np.conj(fra) * frb
    res = np.fft.ifft(sf, axis=0)
    if len(a.shape) > 1:
        # average over all particles/molecules
        cor = (np.real(res[:len(a)]) / np.array(range(len(a), 0, -1))[:, np.newaxis]).mean(axis=1)
    else:
        cor = np.real(res[:len(a)]) / np.array(range(len(a), 0, -1))
    return cor


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
