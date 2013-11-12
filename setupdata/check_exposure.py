#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################
# PROGRAM NAME
#
#
##############################
"""
Created on Fri Nov  1 11:10:44 2013

@author: clackner

program details
"""
import sys, os, re
import argparse
import numpy as np
import glob

import matplotlib.figure as figure
from matplotlib.backends.backend_ps import FigureCanvasPS as FigCanvasPS
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvasA

import pyfits

def main():

    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("infile", help='input file')
    args = parser.parse_args()

    data = pyfits.open(args.infile)[1].data
    files = pyfits.open('/data/bcg_img/redo_oct2013/redo_input.fits')[1].data
    
    for d in data:
        imghd = pyfits.open('../bcg_data/images/'+d['FILENAME']+'.fits')[0].header
    
        if d['IDENT'] in files['IDENT']:
            xx=np.where(files['IDENT']==d['IDENT'])[0][0]
            imghd = pyfits.open('../bcg_data/redo_oct2013/'+files['FILENAME'][xx])[0].header
            #print files[xx]['IDENT'], files[xx]['FILENAME'], imghdf['EXPTIME']
        exptime=imghd['EXPTIME']
        if exptime!=4056.0:
            print d['IDENT'], d['FILENAME'], exptime #imghd['EXPTIME']
    

    return 0


if __name__=='__main__':
    main()
    
