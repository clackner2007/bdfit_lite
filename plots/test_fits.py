#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 14:43:07 2013

@author: clackner
"""

import sys
import pyfits
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

def totalflux(profile, cutoff=True):

    n=profile[2]
    k = n*np.exp(0.6950-0.1789/n)
    reff=profile[1]
    q=profile[3]

    flux = 2.0 * np.pi * reff**2 * profile[0] * np.exp(k) * \
        n * k**(-2.0*n) * gamma(2.*n) * q


    if cutoff:
        if( (n - 1.0) < 1.0e-5 ):
            if( n == 1.0 ):
                flux *= 0.980472
            else:
                flux *= 0.98042454 - 0.06495062*(n-1.) - \
                0.02627548*(n-1.)**2 + 0.02277057*(n-1.)**3 - \
                0.21949023*(n-1)**4 - 0.34741306*(n-1.)**5
        else:
            if( n == 4.0 ):
                flux *= 0.93666
            else:
                flux *= 0.936587 - 2.7233496e-2*(n-4.) + \
                1.4125945e-3*(n-4.)**2 + 4.06329317e-4*(n-4.)**3 - \
                1.82907176e-4*(n-4.)**4 +2.2355412e-5*(n-4.)**5

    return flux


def mag(flux, zeropt=25.959):
    return -2.5*np.log10(flux) + zeropt

def main():

    f1='../bcg_data/imgs_BCG_comb.fits'
    f2='../bcg_data/outputs/fullimg_totalRAW.fits'
    f3='../bcg_data/outputs/trunc_totalRAW.fits'

    dgalfit = pyfits.open(f1)[1].data
    dfull = pyfits.open(f2)[1].data
    dtrunc = pyfits.open(f3)[1].data

    galid = list(np.intersect1d(dfull.IDENT, dtrunc.IDENT))
    idgalfit = np.array([dict(zip(dgalfit.IDENT, range(len(dgalfit))))[i]
                        for i in galid])
    idfull = np.array([dict(zip(dfull.IDENT, range(len(dfull))))[i]
                        for i in galid])
    idtrunc = np.array([dict(zip(dtrunc.IDENT, range(len(dtrunc))))[i]
                        for i in galid])
    dgalfit=dgalfit[idgalfit]
    dfull=dfull[idfull]
    dtrunc=dtrunc[idtrunc]

    print dgalfit.IDENT
    print dfull.DOF_SERSIC
    print dtrunc.DOF_SERSIC

    print dfull['SERSICFIT'][np.where(dfull.IDENT==192620)]
    print dfull['XLEN'][np.where(dfull.IDENT==192620)],
    print dfull['YLEN'][np.where(dfull.IDENT==192620)]
    print dfull['FILENAME'][np.where(dfull.IDENT==192620)]

    print dtrunc['IDENT'][np.where(dtrunc['SERSICFIT'][:,1]-dgalfit['RE_SERSIC'] > 50)]
    print dtrunc['SERSICFIT'][np.where(dtrunc['SERSICFIT'][:,1]-dgalfit['RE_SERSIC'] > 50),2]
    print dfull['SERSICFIT'][np.where(dtrunc['SERSICFIT'][:,1]-dgalfit['RE_SERSIC'] > 50),2]
    print dgalfit['N_SERSIC'][np.where(dtrunc['SERSICFIT'][:,1]-dgalfit['RE_SERSIC'] > 50)]



    plt.ion()

    #dev size compare
    plt.plot(dgalfit['RE_DEV'], dfull['DVCFIT'][:,1]-dgalfit['RE_DEV'], 'r.', label='full')
    plt.plot(dgalfit['RE_DEV'], dtrunc['DVCFIT'][:,1]-dgalfit['RE_DEV'], 'bx', label='trunc')
    plt.xlabel('galfit')
    plt.ylabel('cnl-galfit')
    plt.legend()
    plt.savefig('dev_re_comp.png')
    #
    #ser size compare
    plt.clf()
    plt.subplot(221)
    plt.plot(dgalfit['RE_SERSIC'], dfull['SERSICFIT'][:,1]-dgalfit['RE_SERSIC'], 'r.', label='full')
    plt.plot(dgalfit['RE_SERSIC'], dtrunc['SERSICFIT'][:,1]-dgalfit['RE_SERSIC'], 'bx', label='trunc')
    plt.xlabel('galfit re')
    plt.ylabel('cnl-galfit re')
    plt.legend()
    plt.subplot(222)
    plt.plot(dgalfit['N_SERSIC'], dfull['SERSICFIT'][:,2]-dgalfit['N_SERSIC'], 'r.', label='full')
    plt.plot(dgalfit['N_SERSIC'], dtrunc['SERSICFIT'][:,2]-dgalfit['N_SERSIC'], 'bx', label='trunc')
    plt.xlabel('galfit n')
    plt.ylabel('cnl-galfit n')
    plt.subplot(223)
    plt.plot(dgalfit['Q_SERSIC'], dfull['SERSICFIT'][:,3]-dgalfit['Q_SERSIC'], 'r.', label='full')
    plt.plot(dgalfit['Q_SERSIC'], dtrunc['SERSICFIT'][:,3]-dgalfit['Q_SERSIC'], 'bx', label='trunc')
    plt.xlabel('galfit q')
    plt.ylabel('cnl-galfit q')
    plt.subplot(224)
    plt.plot(dgalfit['MAG_SERSIC'],
             mag([totalflux(dd) for dd in dfull['SERSICFIT']])-
             dgalfit['MAG_SERSIC'], 'r.', label='full')
    plt.plot(dgalfit['MAG_SERSIC'],
             mag([totalflux(dd) for dd in dtrunc['SERSICFIT']])-
             dgalfit['MAG_SERSIC'], 'bx', label='trunc')
    plt.xlabel('galfit mag')
    plt.ylabel('cnl-galfit mag')
    plt.savefig('ser_comp.png')

    #
    #two component fits

    return 0



if __name__=='__main__':
    main()