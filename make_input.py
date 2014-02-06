#!/usr/bin/env python
# -*- coding: utf-8 -*-
###########################
#
# make_input.py
#
# this is an example of how to make input for fit_sample.pro
#
##############################
"""
Created on Thu Feb  6 11:39:36 2014

@author: clackner

program details
"""
import sys, os, re
import argparse
import numpy as np
import pyfits


def main():

    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-t", "--temp",
               action='store_true', default=False, help="temp opt",
               dest='temp')
    args = parser.parse_args()


    #list of names, atlas-ids, parent-ids
    #you would read this in from someplace else
    
    #nsa galaxies use IAU names
    names = ['J011058.90+330907.9']
    atlasids = ['0']
    parentids = ['92']
    
    
    #now, if you don't want to use the default profiles with nothing held fixed
    #you need to set the input parameters here
    #the things to set are NAMEOFPROFILE_FIX, NAMEOFPROFILE_VAL
    #there should be arrays of the size of the output 
    #parameters (8x the number of profiles)
    #The name of the profile is whatever you want, but you need to give it
    #to fit_sample, so it can find the input parameters
    #XXX_FIX is a boolean array which =1 for fixed parameters
    #XXX_VAL is the input values, including fixed values
    #if you don't fix the initial position of flux, default values will 
    #be chosen by the code (based on image size)
    #the parameters (in order) are:
    #0--the surface brightness
    #1 --the half light radius (refers to the semimajor axis)
    #2 --Sersic index
    #3 --axis ratio (minor/major)
    #4 --shape of isophote, default is zero and fixed, 
    #   but if you specify XXX_FIX you need to fix the 4th parameter
    #5 --x coordinate of center of profile, be default multiple profiles 
    #    share a center
    #6 --y coordinate of center
    #7 --position angle (radians counterclockwise from x-axis)
    
    #first profile will be de Vauc with only the isophote shape fixed
    DVC_FIX = np.ndarray((len(names), 8))
    DVC_VAL = np.ndarray((len(names), 8))
    #second profile will be EXP+DVC with all the de Vauc. parameters except
    #the normalization fixed
    EXPDVC_FIX = np.ndarray((len(names), 16))
    EXPDVC_VAL = np.ndarray((len(names), 16))
    for i in range(len(names)):
        DVC_FIX[i,:] = [0,0,1,0,1,0,0,0]
        #you need to pick starting values for the size, and the axis ratio
        #starting values for the surface brightness and the position will
        #be set in the code, don't set the flux to zero, as the code
        #just rescales values (so can't rescale 0)
        DVC_VAL[i,:] = [1.0,10.,4.0,0.7,0.0,0.0,0.0,0.1]
        #for two component fits, always but the 'bulge' (smaller +higher sersic)
        #profile first
        EXPDVC_FIX[i,:] = [0,1,1,1,1,1,1,1,0,0,1,0,1,0,0,0]
        #you'll have to source the fixed values from somewhere, either
        #based on the radius, or based on previously fitting the profiles
        EXPDVC_VAL[i,:] = [1.0, 78.9, 4.0, 0.75, 0.0, 715.217, 725.011, 2.1801,
                            0.0, 100.0, 1.0, 0.75, 0.0, 0.0, 0.0, 2.2]
                            
         #it might be good to have different input scripts, for example
         #you could make 1 input file to fit the sersic/dvc profile and after
         #that runs, make another file for the two component fits, which 
         #takes the olds fits as inputs for the VALs arrays
    


    #put everything into FITS format
    col1 = pyfits.Column(name='NAME', format='A19', array=names)
    col2 = pyfits.Column(name='ATLAS_ID', format='J', array=atlasids)
    col3 = pyfits.Column(name='PARENT_ID', format='J', array=parentids)
    col4 = pyfits.Column(name='DVC_FIX', format='8D', array=DVC_FIX)
    col5 = pyfits.Column(name='DVC_VAL', format='8D', array=DVC_VAL)
    col6 = pyfits.Column(name='EXPDVC_FIX', format='16D', array=EXPDVC_FIX)
    col7 = pyfits.Column(name='EXPDVC_VAL', format='16D', array=EXPDVC_VAL)
    
    cols = pyfits.ColDefs([col1, col2, col3, col4,
                           col5, col6, col7])
    
    hdu = pyfits.new_table(cols)
    hdu.writeto('sample_input.fits', clobber=True)
    
   
    return 0


if __name__=='__main__':
    main()
    
