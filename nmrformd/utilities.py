#!/usr/bin/env python3
"""Utilities for NMRForMD package."""
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2022 Authors and contributors
# Simon Gravelle
#
# Released under the GNU Public Licence, v3 or any higher version
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
from scipy import constants as cst


def fourier_transform(a):
    """
    Calculate Fourier transform.

    Wrap function that takes the data in real space with
    columns time, frequency and returns the data in Fourier space
    with columns frequency, signal
    Units in : ps, signal
    Units out : MHz, s*signal

    Credit to the fourier transform function of MAICoS
    https://maicos-devel.gitlab.io/maicos/index.html
    """
    _dt = (a[1, 0] - a[0, 0]) * cst.pico  # second
    return np.vstack((np.fft.rfftfreq(len(a), _dt)
                     / cst.mega, np.fft.rfft(a[:, 1], norm=None) * _dt * 2)).T

def correlation_function(data):
    """Use fast fourier transforms to calculate the correlation function.

    Credit to the correlation function of MAICoS
    https://maicos-devel.gitlab.io/maicos/index.html
    """
    data_0 = np.append(data, np.zeros(2 ** int(np.ceil((np.log(len(data))
                                   / np.log(2)))) - len(data)), axis=0)
    data_1 = np.append(data_0, np.zeros(data_0.shape), axis=0)
    fra = np.fft.fft(data_1, axis=0)
    sf = np.conj(fra) * fra
    res = np.fft.ifft(sf, axis=0)
    return np.real(res[:len(data)]) / np.array(range(len(data), 0, -1))

def find_nearest(array, value):
    """Find nearest value."""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
