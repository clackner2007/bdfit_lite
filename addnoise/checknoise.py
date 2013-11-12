#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################
# checknoise.py
#
##############################
"""
Created on Tue Nov 12 13:28:51 2013

@author: clackner

program details
"""
import sys, os, re
import argparse
import numpy as np
import scipy.special

import pyfits

import matplotlib.figure as figure
from matplotlib.backends.backend_ps import FigureCanvasPS as FigCanvasPS
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvasA


def totalFlux(data, profile='DVC'):
    """
    get the total profile flux
    """
    prof = data[profile+'FIT']
    k = prof[:,2]*np.exp(0.695 - 0.1789/prof[:,2])
    flux1 = 2*np.pi*prof[:,1]**2 * prof[:,0] * np.exp(k) * prof[:,2] * \
        k**(-2.0*prof[:,2]) * scipy.special.gamma(2.0*prof[:,2])*prof[:,3]
        
    if profile in ['DVC', 'SERSIC']:
        return flux1
    else:
        k = prof[:,10]*np.exp(0.695 - 0.1789/prof[:,10])
        flux2 = 2*np.pi*prof[:,9]**2 * prof[:,8] * np.exp(k) * prof[:,10] * \
        k**(-2.0*prof[:,10]) * scipy.special.gamma(2.0*prof[:,10])*prof[:,11]
        return flux1 + flux2

def getErr(data, profile='DVC'):
    """
    print out error information for profile
    assume all entris in data refer to one galaxy
    """
    
    print "mean size", np.mean(data[profile+'FIT'][:,1]), np.std(data[profile+'FIT'][:,1]),
    print np.mean(data['PERR_'+profile][:,1])
    if profile in ['SERSIC', 'DSERSIC']:
        print "mean n", np.mean(data[profile+'FIT'][:,2]), np.std(data[profile+'FIT'][:,2]),
        print np.mean(data['PERR_'+profile][:,2])
    if profile in ['DDVC', 'DSERSIC']:
        print "mean size", np.mean(data[profile+'FIT'][:,9]), 
        print np.std(data[profile+'FIT'][:,9]),
        print np.mean(data['PERR_'+profile][:,9])
        if profile in ['DSERSIC']:
            print "mean n", np.mean(data[profile+'FIT'][:,10]), 
            print np.std(data[profile+'FIT'][:,10]),
            print np.mean(data['PERR_'+profile][:,10])
        print "mean B/T", np.mean(data['FLUX_RATIO_'+profile]),
        print np.std(data['FLUX_RATIO_'+profile])
        print "mean reff tot", np.mean(data['REFF_'+profile]), np.std(data['REFF_'+profile])
    print 'mean flux', np.mean(totalFlux(data, profile)), np.std(totalFlux(data, profile))
    

def main():

    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("imlist", help='image list')
    parser.add_argument("-e", "--eps",
                 action='store_true', default=False, help='make eps plots',
                 dest='epsPlot')
    args = parser.parse_args()
    
    data = pyfits.open(args.imlist)[1].data
    ids = set(data['IDENT'])
    
    for ident in ids:
        keep = np.where(data['IDENT']==ident)[0]
        print ident#, data[keep]['DDVCFIT'][:,[1,9]]
        keep = keep[1:]
        print 'DVC'
        getErr(data[keep], 'DVC')
        print 'SERSIC'
        getErr(data[keep], 'SERSIC')
        print 'DSERSIC'
        getErr(data[keep], 'DSERSIC')
    

    FigCanvas = FigCanvasPS if args.epsPlot else FigCanvasA
    ending='.eps' if args.epsPlot else '.png'

    fig=figure.Figure()
    canv=FigCanvas(fig)

    return 0


if __name__=='__main__':
    main()
    
