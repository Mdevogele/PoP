#!/usr/bin/env python

""" PoP_Toolbox - Function toolbox for the Polarimetry pipeline
    v1.0: 2018-04-17, mdevogele@lowell.edu
"""

import numpy as np

from astropy.io import fits

from PoP_CheckInstrument import CheckInstrument


def GetData(Filename):
    telescope, obsparam = CheckInstrument(Filename)
    hdulist = fits.open(Filename)
    if 'NOT' in telescope:
        data = hdulist[1].data
    else:
        data = hdulist[0].data
        
    return data

def SetData(Filename,data):
    telescope, obsparam = CheckInstrument(Filename)
    hdulist = fits.open(Filename)
    if 'NOT' in telescope:
        hdulist[1].data = data
    else:
        hdulist[0].data = data
        
    return hdulist
        
