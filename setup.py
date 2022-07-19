#!/usr/bin/env python3
"""Setup file for NMRforMD package."""
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2022 Authors and contributors
# Simon Gravelle
#
# Released under the GNU Public Licence, v3 or any higher version
# SPDX-License-Identifier: GPL-3.0-or-later

from setuptools import setup

setup(name='nmrformd',
      version='v0.0.9',
      description='Calculate NMR relaxation time from \
                   molecular dynamics trajectory file',
      long_description=open('README.rst').read(),
      long_description_content_type='text/x-rst',
      url='https://github.com/simongravelle/nmrformd',
      download_url='https://github.com/simongravelle/nmrformd/archive/refs/tags/v0.0.9.tar.gz',  # noqa
      author='Simon Gravelle',
      author_email='simon.gravelle@live.fr',
      license='GNU GENERAL PUBLIC LICENSE',
      packages=['nmrformd'],
      zip_safe=False,
      install_requires=[
       "mdanalysis",
       "pytest",
       "numpy",
       "scipy",
      ]
      )
