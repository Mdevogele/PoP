#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
_PoP_conf configuration file for Polarimetry Pipeline

2019-04-01, mdevogele@lowell.edu
"""
from __future__ import print_function

import os
import sys
import warnings

try:
#    from astropy import wcs
    from astropy.io import fits
except ImportError:
    print('Module astropy not found. Please install with: pip install astropy')
    sys.exit()

try:
    import numpy as np
except ImportError:
    print('Module numpy not found. Please install with: pip install numpy')
    sys.exit()
    

# potential FITS header keywords for looking up the instrument
# any unique header keyword works as a potential identifier
instrument_keys = ['INSTRUME', 'LCAMMOD']

# suppress runtime and astropy warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

# only import if Python3 is used
if sys.version_info > (3, 0):
    from past.builtins import execfile


# read polarimetry pipeline root path from environment variable
rootpath = os.environ.get('POPPIPEDIR')
if rootpath is None:
    print('ERROR: POPPIPEDIR variable has not been set')
    sys.exit(0)


execfile(rootpath+'/setup/telescopes.py')
