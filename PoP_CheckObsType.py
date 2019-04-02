#!/usr/bin/env

""" SP_ObsTypeCheck - Script that check the observation type for each file and correct it if necessary
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

import _PoP_conf
import sys
import logging
import numpy as np
import os
import shutil

import re

from astropy.io import fits

from PoP_CheckInstrument import CheckInstrument

from astroquery.jplhorizons import Horizons
from astroquery.simbad import Simbad

import astropy.coordinates as coord
from astropy import units as u

_PoP_conf.filenames = []


def CheckObsType(filenames):
    
    List_ID_LoPol = ['GD 319', 'BD+33  2642']
    
    List_ID_LoPol_Comp = {'GD  319': 'GD319', 'BD+33  2642': 'BD332642'}    
    
    List_ID_HiPol = ['HD 155197', 'HD 155528']
    
    List_ID_HiPol_Comp = {'HD 155197': 'HD155197', 'HD 155528': 'HD155528'}    
    
    telescope, obsparam = CheckInstrument(filenames)
    
    for idx, filename in enumerate(filenames):
        try:
            hdulist = fits.open(filename, mode='update' ,ignore_missing_end=True)
        except IOError:
            logging.error('cannot open file %s' % filename)
            print('ERROR: cannot open file %s' % filename)
            filenames.pop(idx)
            continue

        if 'NOT' in telescope:
            data = hdulist[1].data
        else:
            data = hdulist[0].data
            
        Med = np.nanmedian(data[10:-10,10:-10])
        Max = np.nanmax(data[10:-10,10:-10])
        std0 = np.nanstd(np.median(data[10:-10,10:-10],axis=0))
        std1 = np.nanstd(np.median(data[10:-10,10:-10],axis=1))
        
        if telescope == '2.12m':
            TIME = hdulist[0].header[obsparam['date_keyword']].replace('-','').replace(':','').split('.')[0]
            EXPTIME = str(hdulist[0].header[obsparam['exptime']]).replace('.','s')
            FILTER = str(hdulist[0].header[obsparam['filter']])
            if 'image' in hdulist[0].header[obsparam['obstype']]:
                TYPE = 'Image'
                Name= telescope + '_' + 'POLIMA2' + '_' + FILTER + '_' + TIME + '_Image_' + EXPTIME + '.fits'
                hdulist.close()
                _PoP_conf.filenames.append(Name)
                print('%s changed to %s' % (filename, Name))
                shutil.copy(filename,Name)
            if 'flat' in hdulist[0].header[obsparam['obstype']]: # Acquisition files
                print('Flat image')
                TYPE = 'FLAT'
                Name= telescope + '_' + 'POLIMA2' + '_' + FILTER + '_' + TIME + '_FLAT_' + EXPTIME + '.fits'
                hdulist.close()
                shutil.copy(filename,Name)
                _PoP_conf.filenames.append(Name)
                print('%s changed to %s' % (filename, Name))
            if 'zero' in hdulist[0].header[obsparam['obstype']]:
                print('BIAS image')
                TYPE = 'BIAS'
                Name= telescope + '_' + 'POLIMA2' + '_' + TIME + '_BIAS_' + EXPTIME + '.fits'
                hdulist.close()
                shutil.copy(filename,Name)
                _PoP_conf.filenames.append(Name)
                print('%s changed to %s' % (filename, Name))
                
            if 'object' in hdulist[0].header[obsparam['obstype']]:
                
                RotAng = np.abs(float(hdulist[0].header[obsparam['polang']]))
                TYPE = 'Unknown'
                hdulist[0].header[obsparam['obstype']] = 'OBJECT'
                OBJECT = hdulist[0].header[obsparam['object']].replace(' ','').replace('/','')
                print(obsparam['object'])
                
                SEP = re.findall('\d+|\D+',OBJECT.replace('(','').replace(')',''))[0]
            
                if OBJECT.replace(' ',''):
                    print(OBJECT)
                    try:
                        Horizons(id=OBJECT).ephemerides()
                        TYPE = 'Asteroid'
                    except ValueError:
                        pass    
                    try:
                        Horizons(id=SEP).ephemerides()
                        TYPE = 'Asteroid'
                    except ValueError:
                        pass
                    
                    if TYPE == 'Unknown':
                        try: 
                            #result_table = Simbad.query_object(OBJECT)
                            result_table = Simbad.query_region(coord.SkyCoord(hdulist[0].header[obsparam['ra']] + ' ' + hdulist[0].header[obsparam['dec']] ,unit=(u.hourangle,u.deg), frame='icrs'), radius='0d30m00s')
                            T = result_table['MAIN_ID']
                            for elem in T:
                                if elem in List_ID_LoPol: 
                                    OBJECT = List_ID_LoPol_Comp[elem]
                                    TYPE = 'LoPol'
                            for elem in T:
                                if elem in List_ID_HiPol: 
                                    OBJECT = List_ID_HiPol_Comp[elem]
                                    TYPE = 'HiPol'        
                        except:
                            pass
                        if TYPE == 'Unknown':
        #                       except TypeError:
                            if '(' in OBJECT and ')' in OBJECT:
                                Number = OBJECT[OBJECT.find("(")+1:OBJECT.find(")")]
                                if Number.isdigit():
                                    try:
                                        Horizons(id=Number).ephemerides()
                                        TYPE = 'Asteroid'
                                    except ValueError:
                                        TYPE = 'Unknown'
                                    if OBJECT[0:3].isdigit():
                                        if not ' ' in OBJECT:
                                            Name = OBJECT[0:4] + ' ' + OBJECT[4:]
                                            try:
                                                Horizons(id=Name).ephemerides()
                                                TYPE = 'Asteroid'
                                            except ValueError:
                                                TYPE = 'Unknown'
                                                            
                                                        
                    Name= telescope  + '_POLIMA2_' + TIME + '_' + FILTER + '_' + TYPE + '_' + OBJECT + '_' + str(RotAng) + '_' + EXPTIME + '.fits'
                    
                    shutil.copy(filename,Name)
                                                    
                    _PoP_conf.filenames.append(Name)
                    print('%s changed to %s' % (filename, Name))
                    hdulist.close()    
                       

                
    for elem in filenames:
        os.remove(elem)                  
       
    return None
    
    