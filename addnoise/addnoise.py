#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################
# addnoise.py
#
##############################
"""
Created on Tue Nov 12 08:52:02 2013

@author: clackner

program details
"""
import sys, os, re, copy
import argparse
import numpy as np

import pyfits

import matplotlib.figure as figure
from matplotlib.backends.backend_ps import FigureCanvasPS as FigCanvasPS
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvasA


def addnoise(image, ivar, scaleup=2.0):
    """
    add noise by a factor of scaleup
    """
    
    sigma = np.sqrt(1.0/np.median(ivar))
    new_ivar = 1.0/(scaleup**2) * ivar
    new_image = copy.copy(image)

    noise = np.random.normal(0, np.sqrt(scaleup**2-1.0)*sigma, image.shape)
    new_image += noise

    return (new_image, new_ivar)    


def getImage(ident, filename, path='/data/bcg_img/', header=True):
    
    img = pyfits.open(path+'images/'+filename+".fits")[0].data
    ivar = pyfits.open(path+'ivar/'+filename+'.wht.fits')[0].data
    imhead = pyfits.open(path+'images/'+filename+".fits")[0].header
    ivhead = pyfits.open(path+'ivar/'+filename+'.wht.fits')[0].header
    if not header:
        return img, ivar
    else:
        return img, ivar, imhead, ivhead

def main():

    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("imlist", help='image list')
    parser.add_argument('nimg', help='number of new images', type=int)
    parser.add_argument('-s', '--sigma', help='multiplier for noise',
                        dest='scale', default=2.0, type=float)
    parser.add_argument('-o', '--outpath', dest='outpath',
                        default='', help='output image path')
    parser.add_argument('-i', '--inpath', dest='inpath',
                        default='/data/bcg_img/')
    parser.add_argument("-e", "--eps",
                 action='store_true', default=False, help='make eps plots',
                 dest='epsPlot')
    args = parser.parse_args()
    
    data = pyfits.open(args.imlist)[1].data
    ids = []
    files = []
    for datum in data[144:145]:
        img, ivar, imhead, ivhead = getImage(datum['IDENT'], datum['FILENAME'],
                                             path=args.inpath)
        ids.append(datum['IDENT'])
        files.append(datum['FILENAME'])
        for i in range(args.nimg):
            img2, ivar2 = addnoise(img, ivar, scaleup=args.scale)
            hdu1=pyfits.PrimaryHDU(img2, header=imhead)
            hdu1.writeto(args.outpath+'images/'+
                        datum['FILENAME']+'_%02d.fits'%i, clobber=True)
            hdu2=pyfits.PrimaryHDU(ivar2, header=ivhead)
            hdu2.writeto(args.outpath+'ivar/'+
                        datum['FILENAME']+'_%02d.wht.fits'%i, clobber=True)
            ids.append(datum['IDENT'])
            files.append(datum['FILENAME']+'_%02d'%i)

    hdu = pyfits.new_table((pyfits.Column(name='IDENT',format='J',
                                          array=np.asarray(ids)),
                            pyfits.Column(name='FILENAME',format='40A',
                                          array=np.asarray(files))))
    hdu.writeto(args.outpath+'add_noise_list.fits')
    return 0


if __name__=='__main__':
    main()
    
