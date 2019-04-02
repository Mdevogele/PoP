#!/usr/bin/env python

""" PoP_Flat - wrapper for creating flat(s) 
    v1.0: 2018-04-17, mdevogele@lowell.edu
"""

import argparse, shlex

import PoP_toolbox as tb
import numpy as np

from astropy.io import fits

from itertools import groupby
from operator import itemgetter

import PoP_diagnostics as diag
   

def rebin(arr, new_shape):
    """Rebin 2D array arr to shape new_shape by averaging."""
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)


class Result(object):
    def Indiv_to_Master(self, method, Flat):
        Master = getattr(self, method, lambda x: getattr(self, 'Default')(Flat))(Flat)
        return Master
 
    def Median(self,Flat):
        FF = []
        for elem in Flat:
            FF.append(elem/np.median(elem))
            Master = np.median(FF,axis=0)/np.median(FF)
        return Master
 
    def Mean(self,Flat):
        Master = np.mean(Flat,axis=0)/np.mean(Flat)
        return Master

    def Default(self, Flat):
        Master = np.median(Flat,axis=0)/np.median(Flat)
        print("Invalid method, use of the median as default")
        return Master 




def Create_Flat(filenames,MasterName,Verbose,Bias,Diagnostic,Method):


    if Bias:
        if Verbose: 
            print('Opening master bias file')            
        Bias_data = tb.GetData(Bias)
    
    
    if Verbose:
        print('Beginning flat processing')
        print('Processing files:')
        print('index \t filename')
        for idx,elem in enumerate(filenames):
            print('{} \t {}'.format(idx+1,elem))
        print('Is using {} as master bias'.format(Bias))
        
                
    Flat = []
    for idx,elem in enumerate(filenames):
        print('{} \t {}'.format(idx+1,elem))
        data = tb.GetData(elem)
        Flat.append(data - Bias_data)

    print('Creating the master flat')

    Res = Result()
    MasterFlat = Res.Indiv_to_Master(Method, Flat)   

    MasterFlat[MasterFlat<0.8] = 1
    MasterFlat[MasterFlat>1.2] = 1

    hdulist = tb.SetData(Bias,np.float32(MasterFlat))

    hdulist[0].header['BZERO'] = 0
    hdulist.writeto(MasterName, overwrite = True)
    hdulist.close()    


    if Diagnostic:
        diag.create_website('Flat_Log.html')
        diag.add_BiasSummary(filenames,MasterName,'Flat_Log.html')
        diag.add_FlatList(filenames,'Flat_Log.html')

    if Verbose:
        print('Master flat save to {}'.format(MasterName))
        hdulist = fits.open(MasterName)
        data = hdulist[0].data
        print('Statistics of the master flat')
        print('Mean: {} \t Median: {} \t std: {}'.format(np.mean(data), np.median(data), np.std(data)))
        print('End of flat processing')
        hdulist.close()
            

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Processing and creation of master flats')


    parser.add_argument('-b',
                        help='Name of the master bias to use \n || Can use None if no bias are available',
                        default='MasterBias.fits')
    parser.add_argument('-v',
                        help='Increase verbosity',
                        action="store_true")    
    parser.add_argument('-o',
                        help='Prefix of the master flats files',
                        default='MasterFlat.fits')  
    parser.add_argument('-m',
                        help='Method to use to compute the master bias: Mean, Median',
                        default='Median') 
    parser.add_argument('-d',
                        help='Enable or disable the diagnostic',
                        action="store_false",
                        default = True)  
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    args = parser.parse_args()
    Bias = args.b
    Verbose = args.v
    MasterName = args.o
    Method = args.m 
    filenames = args.images    
    Diagnostic = args.d

    
    # call run_the_pipeline only on filenames
    Create_Flat(filenames,MasterName,Verbose,Bias,Diagnostic,Method)
    pass
