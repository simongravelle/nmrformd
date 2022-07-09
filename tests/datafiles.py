#!/usr/bin/env python
"""Import datafiles."""
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2022 Authors and contributors
# Simon Gravelle
#
# Released under the GNU Public Licence, v3 or any higher version
# SPDX-License-Identifier: GPL-3.0-or-later

from pkg_resources import resource_filename

SPCE_ITP = resource_filename(__name__, '../data/spce.itp')
SPCE_GRO = resource_filename(__name__, '../data/spce.gro')