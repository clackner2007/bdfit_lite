#!/usr/bin/env python

##############################
# make-mask.py
#
# Claire Lackner
# 2013-08-16
#
# make masks from sextractor segmentation map
##############################
"""program details"""

import sys, os, re
import argparse
import numpy as np
import pyfits
import copy
import scipy.ndimage as ndimage

def embiggenMask(segmap, objnum, filtersize=10):
    """
    make the masked regions larger
    """
#    import pdb
#    import matplotlib.pyplot as plt
#    pdb.set_trace()
    
    origmask = copy.copy(segmap)
    origmask[segmap==objnum]=0.0
    
    filteredmask = ndimage.filters.maximum_filter(origmask, 
                                                  size=[filtersize, 
                                                        filtersize], 
                                                mode='nearest')
    
    segmap[(filteredmask!=0)] = filteredmask[filteredmask!=0]



def makemask(segfile, bigger=0, clobber=True, path=''):
    #from Kevin's make_mask.pro, just pick the central pixel value as the
    #main object
    seg = pyfits.open(segfile)[0].data
    objmask = seg[seg.shape[0]/2, seg.shape[1]/2]
    if objmask==0:
        print 'center pixel has no object'
        return 0    
    
    #make mask larger
    if bigger:
        embiggenMask(seg, objmask, filtersize=bigger)
    
    mask = seg
    mask[mask==objmask] = 0
    mask[mask!=0] = 1

    maskHDU = pyfits.PrimaryHDU(mask)
    maskHDU.writeto(path+os.path.basename(segfile.split('_')[-2]+'_mask.fits'),
                    clobber=clobber)

def main():

    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('segfile', help='segmentation map')
    parser.add_argument('-b', '--big', help='make mask larger',
                        default=0, type=int, dest='bigger')
    parser.add_argument('-c', '--clobber', action='store_true',
                        default=False, dest='clobber')
    args = parser.parse_args()
    makemask(args.segfile, args.bigger, args.clobber)
    

    return 0


if __name__=='__main__':
    main()
