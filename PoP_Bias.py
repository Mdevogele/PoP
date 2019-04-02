#!/usr/bin/env python

""" PoP_Bias - wrapper for creating bias 
    v1.0: 2018-04-17, mdevogele@lowell.edu
"""

import argparse, shlex

import numpy as np

from astropy.io import fits
from PoP_CheckInstrument import CheckInstrument
import PoP_toolbox as tb

import PoP_diagnostics as diag

def Create_Bias(filenames,MasterName,Verbose,Diagnostic):
    
    print(filenames)
    
    if Verbose:
        print('Beginning bias processing')
        print('Processing files:')
        print('index \t filename')
        for idx,elem in enumerate(filenames):
            print('{} \t {}'.format(idx+1,elem))

  
    if Verbose:
        print('Creating the master bias')
           
    Bias = []
    for image in filenames:
        
        toopen = image
        print(toopen)
        hdulist = fits.open(toopen)
        data = tb.GetData(toopen)
        Bias.append(data)
        
    MasterBias = np.median(Bias,axis = 0 )
    
    hdulist = tb.SetData(image,data)

    hdulist.writeto(MasterName, overwrite = True)
    hdulist.close()
    
    if Verbose:
        print('Master bias save to {}'.format(MasterName))
        hdulist = fits.open(MasterName)
        data = hdulist[0].data
        print('Statistics of the Master bias: ')
        print('Mean: {} \t Median: {} \t std: {}'.format(np.nanmean(data), np.nanmedian(data), np.nanstd(data)))
        print('End of bias processing')

    if Diagnostic:
        diag.create_website('Bias_Log.html')
        diag.add_BiasSummary(filenames,MasterName,'Bias_Log.html')
        diag.add_BiasList(filenames,'Bias_Log.html')        
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Processing and creation of master bias')

    parser.add_argument('-lv',
                        help='decrease verbosity',
                        action="store_false",
                        default = True)    
    parser.add_argument('-o',
                        help='Name of the master bias file',
                        default='MasterBias.fits')  
    parser.add_argument('-d',
                        help='Enable or disable the diagnostic',
                        action="store_false",
                        default = True)   
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()

    Verbose = args.lv
    MasterName = args.o
    filenames = args.images    
    Diagnostic = args.d

    Create_Bias(filenames,MasterName,Verbose,Diagnostic)
    pass