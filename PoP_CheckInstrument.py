#!/usr/bin/env

""" PoP_CheckInstrument - Script that check the instrument for each file
    v1.0: 2019-04-01, mdevogele@lowell.edu
"""
from __future__ import print_function

# Spectroscopic Pipeline
# Copyright (C) 2018  Maxime Devogele, mdevogele@lowell.edu

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.


import sys

import _PoP_conf

from astropy.io import fits


def CheckInstrument(filenames):
    instruments = []
    if isinstance(filenames, str):
        filenames = [filenames]
    
    for idx, filename in enumerate(filenames):
        try:
            hdulist = fits.open(filename, ignore_missing_end=True)
        except IOError:
            print('ERROR: cannot open file %s' % filename)
            continue

        header = hdulist[0].header
        for key in _PoP_conf.instrument_keys:
            if key in header:
                instruments.append(header[key])
                break

    if len(instruments) == 0:
        raise KeyError('cannot identify telescope/instrument; please update' + \
                       '_PoP_conf.instrument_keys accordingly')


    telescope = _PoP_conf.instrument_identifiers[instruments[0]]
    obsparam = _PoP_conf.telescope_parameters[telescope]    
               
    
    
    return telescope, obsparam