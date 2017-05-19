#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2016-2017, Du-lab
# Author: Owen Myers
# Contact: omyers2@uncc.edu
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.

import scipy
import pylab as pl

def gaussian(x, a, b, c):
    return a * pl.exp(-(x-b) ** 2 / (2 * c))
def asym_gaussian(x,a,b,c,d):
    return 2.0*gaussian(x,a,b,c)*0.5*(1+scipy.special.erf(d*(x-b)/pl.sqrt(2)))