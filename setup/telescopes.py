"""
Polarimetric Pipeline Configuation File
2019-04-01, mdevogele@lowell.edu
"""

# Spectroscopic Pipiline
# Copyright (C) 2019  Maxime Devogele, mdevogele@lowell.edu

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

##### telescope/instrument configurations


# Polalima@@San Pedro Martir Mexico
SPM_2m_param = {
    'telescope_instrument' : 'Polima2', # telescope/instrument name
    'telescope_keyword'    : '2.12m',      # telescope/instrument keyword
    'secpix'               : 0.34, # pixel size (arcsec)
                                            # before binning

    # image orientation preferences
    'flipx'                : True,      # Is the wavelength increase with increasing X values ?

    # instrument-specific FITS header keywords
    'binning'              : ('CCDXBIN', 'CCDYBIN'),
                           # binning in x/y, '_blankN' denotes that both axes
                           # are listed in one keyword, sep. by blanks
    'extent'               : ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra'                   : 'RA',  # telescope pointing, RA
    'dec'                  : 'DEC', # telescope pointin, Dec
    'radec_separator'      : ':',   # RA/Dec hms separator, use 'XXX'
                                    # if already in degrees
    'date_keyword'         : 'UTMIDDLE', # obs date/time
                                                  # keyword; use
                                                  # 'date|time' if
                                                  # separate
    'obsmidtime_jd'        : 'MIDTIMJD', # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object'               : 'OBJECT',  # object name keyword
    'filter'               : 'FILTER',  # filter keyword
    'filter_translations'  : {},
                             # filtername translation dictionary
    'exptime'              : 'EXPTIME', # exposure time keyword (s)
    'airmass'              : 'AIRMASS', # airmass keyword
    
    'grating'              : 'GRATING', # grating used
    'obstype'              : 'IMGTYPE', # get the type of observation (bias, flat, object)
    
    'polang'                : 'POLANGLE',
    
    # Limits for BIAS, FLAT, ARCS, OBJECT auto detection
    
    'flat_median'          : 10000,   # if the median of all pixels is > 10000 files considered as a flat
    'bias_std'             : 10,      # if the sum of the std of the median of each axis > 10 files considered as Bias   

}


##### access functions for telescope configurations


implemented_telescopes = ['SPM_2m']

# translate INSTRUME (or others, see _pp_conf.py) header keyword into
# PP telescope keyword
instrument_identifiers = {'Polima2':        '2.12m'
}

# translate telescope keyword into parameter set defined here
telescope_parameters = {'2.12m' :       SPM_2m_param
}


try:
    execfile(rootpath+'/setup/mytelescopes.py')
except IOError:
    pass
