#!/usr/bin/env python

""" PoP_Preproc - Apply the flat and bias to science data
    v1.0: 2019-04-02, mdevogele@lowell.edu
        
"""
import argparse, shlex

import PoP_toolbox as tb

import PoP_diagnostics as diag

import numpy as np


def Preproc(filenames,MasterBias,MasterFlat,Verbose,Suffix,Diagnostic):

    
    Object_out = []
    Object_in = []
    compteur = 0
    for elem in filenames:
        compteur += 1
        Object_out.append(elem.replace('.fits','') + '_' + Suffix)
        Object_in.append(elem[:-5])

    if MasterBias:
        if Verbose: 
            print('Opening master bias file: ' + str(MasterBias))            
        Bias_data = tb.GetData(MasterBias)
        
    if MasterFlat: 
        if Verbose: 
            print('Opening master flat file: ' + str(MasterFlat))            
        Flat_data = tb.GetData(MasterFlat)
    
    for idx, elem in enumerate(filenames):
        data = tb.GetData(elem)
        
        if MasterBias and not MasterFlat:
            data = data.astype(float)- Bias_data.astype(float)
        if MasterBias and MasterFlat:
            data = (data.astype(float) - Bias_data.astype(float))/Flat_data.astype(float)
        if MasterFlat and not MasterBias:
            data = data.astype(float)/MasterFlat.astype(float)
        
        hdulist = tb.SetData(elem,data)
        hdulist[0].header['BZERO'] = 0
        hdulist.writeto(elem.split('.fit')[0] + '_' +  Suffix + '.fits',overwrite = True)
                             
    if Diagnostic: 
        diag.create_website('Pre-processing_Log.html')
        diag.add_PreProc(Object_out,MasterBias,MasterFlat,'Pre-processing_Log.html')
        



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Processing and creation of master flats')
#    parser.add_argument('-prefix', help='data prefix',
#                        default=None)
    parser.add_argument('-s',
                        help='Suffix to add to processed files',
                        default='Proc')
    parser.add_argument('-b',
                        help='Name of the master bias to use \n || Can use None if no bias are available',
                        default='MasterBias.fits')
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_false")    
    parser.add_argument('-f',
                        help='Name of the master flat to use \n || Can use None if no bias are available',
                        default='MasterFlat.fits')    
    
    parser.add_argument('-d',
                        help='Enable or disable the diagnostic',
                        action="store_false",
                        default = True)  
    
    
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
    Suffix = args.s
    MasterBias = args.b
    Verbose = args.v
    MasterFlat = args.f
    filenames = args.images    
    Diagnostic = args.d
    
    print(filenames)
    
    Preproc(filenames,MasterBias,MasterFlat,Verbose,Suffix,Diagnostic)
    pass