#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################
# PROGRAM NAME
#
#
##############################
"""
Created on Fri Nov  1 14:34:04 2013

@author: clackner

program details
"""
import sys, os, re
import argparse
import numpy as np
import pyfits
import glob
from makemask import makemask

def main():

    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('path', help='path to output')
    args = parser.parse_args()
    
    files=glob.glob('../bcg_data/redo_oct2013/images/*3.fits')
    files = sorted(files)
    
    idents=np.ndarray((len(files),), dtype=int)
    filenames=np.asarray([os.path.splitext(os.path.basename(f))[0] 
                            for f in files])
    
    for indfile, f in enumerate(files):
        idents[indfile] = int(os.path.basename(f).split('_')[0])
        print "sex %s -c config.sex -CHECKIMAGE_NAME %s%d_seg.fits"%\
                (f, args.path, idents[indfile])
        os.system("sex %s -c config.sex -CHECKIMAGE_NAME %s%d_seg.fits"%\
                (f, args.path, idents[indfile]))
        makemask('%s%d_seg.fits'%(args.path, idents[indfile]), 0, True, args.path)

    cols=pyfits.ColDefs([pyfits.Column(name='IDENT', format='J',
                                       array=idents),
                        pyfits.Column(name='FILENAME', format='40A',
                                      array=filenames)])
    hdu = pyfits.new_table(cols)
    hdu.writeto(args.path+'redo_input.fits', clobber=True)

    return 0


if __name__=='__main__':
    main()
    
