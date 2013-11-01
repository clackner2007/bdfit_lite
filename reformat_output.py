#!/usr/bin/env python

##############################
# reformat_output.py
#
# Claire Lackner
# 2013-05-14
#
# 
##############################
"""program details"""

import sys, os, re
import argparse
import numpy as np
import pyfits
import scipy.special

def makeCols(names, types, length):
    zeros={'J':int, 'D':float}
    cols=[]
    for i in range(len(names)):
        cols.append(pyfits.Column(name=names[i],
                                  format=types[i],
                                  array=np.zeros(length, 
                                                 dtype=zeros[types[i]])))
    return pyfits.new_table(cols)
                  

def totalflux(a0, r0, n, q, cutoff=True):
    k = n*np.exp(0.6950-0.1789/n)
    flux = 2 * np.pi * r0**2 * a0 * np.exp(k) * n * \
        k**(-2.0*n) * q * scipy.special.gamma(2.0*n)
    if not cutoff:
        return flux
    else:
        factor = np.ones_like(flux)
        print len(flux)
        small = (n-1 < 1.e-5)
        print len(small)
        ns = n[small]
        factor[small] = 0.98042454 - 0.06495062*(ns-1.) - \
            0.02627548*(ns-1.)**2 + 0.02277057*(ns-1.)**3 - \
            0.21949023*(ns-1)**4 - 0.34741306*(ns-1.)**5 
        
        large = ~small
        nl = n[large]
        factor[large] = 0.936587 - 2.7233496e-2*(nl-4.) + \
            1.4125945e-3*(nl-4.)**2 + 4.06329317e-4*(nl-4.)**3 - \
            1.82907176e-4*(nl-4.)**4 +2.2355412e-5*(nl-4.)**5
        if (np.abs(n-1) < 1.0e-5).any():
            factor[np.abs(n-1) < 1.0e-5] = 0.980472
        if (np.abs(n-4) < 1.0e-5).any():
            factor[np.abs(n-4) < 1.0e-5] = 0.933666        
        return flux*factor


def main():

    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("infile", help='input file')
    parser.add_argument('path', help='path to output')
    parser.add_argument('-n', '--new', dest='new', help='new output files',
                        action='store_true', default=False)
    args = parser.parse_args()
    
    data = pyfits.open(args.infile)[1].data

    #make output files for each type of fit
    outcols = ['IDENT', 'R_INNER', 'N_INNER', 'R_OUTER', 'N_OUTER',
               'INNER_FLUX_FRAC', 'R_TOT', 'STATUS', 'Q_INNER',
               'Q_OUTER', 'Q_TOT', 'R_INNER_ERR', 'R_OUTER_ERR', 
               'R_TOT_ERR', 'TOTAL_FLUX', 'MAD_RATIO', 'RED_CHI2',
               'I0_INNER', 'I0_OUTER']
    outtype = ['J', 'D', 'D', 'D', 'D', 'D', 'D', 'J', 'D', 'D', 'D',
               'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D']
    

    if args.new:
        outfiles = ['dvcsersic', 'ddvc', 'dvcexp']
        innames = ['DVCSERSIC', 'DDVC', 'DVCEXP']
        innums = [0, 1, 2]
    else:
        outfiles = ['sersic', 'dblsersic', 'dbldvc', 'expsersic', 'dvc']
        innames = ['SERSIC', 'DSERSIC', 'DDVC', 'EXPSERSIC', 'DVC']
        innums = [1, 0, 4, 3, 2]
    
    nentry = len(data)
    
    for n, name in enumerate(innames):
        newtab = makeCols(outcols, outtype, nentry)
        #print len(data['IDENT']), len(newtab.data['IDENT'])
        newtab.data.field('IDENT')[:] = data['IDENT']

        if (name=='SERSIC') or (name=='DVC'):
            newtab.data.field('R_INNER')[:] = data[name+'FIT'][:,1]
            newtab.data.field('R_INNER_ERR')[:] = data['PERR_'+name][:,1]
            newtab.data.field('Q_INNER')[:] = data[name+'FIT'][:,3]
            newtab.data.field('N_INNER')[:] = data[name+'FIT'][:,2]
            newtab.data.field('INNER_FLUX_FRAC')[:] = 1.0
            newtab.data.field('R_TOT')[:] = newtab.data.field('R_INNER')
            newtab.data.field('R_TOT_ERR')[:] = data['PERR_'+name][:,1]
            newtab.data.field('I0_INNER')[:] = (data[name+'FIT'][:,0]*
                                                np.exp(data[name+'FIT'][:,2]*
                                                       np.exp(0.6950-0.1789/
                                                              data[name+'FIT'][:,2])))
            newtab.data.field('TOTAL_FLUX')[:] \
                = totalflux(data[name+'FIT'][:,0], 
                            data[name+'FIT'][:,1],
                            data[name+'FIT'][:,2],
                            data[name+'FIT'][:,3])

        else:
            newtab.data.field('R_INNER')[:] = data[name+'FIT'][:,1+8]
            newtab.data.field('R_INNER_ERR')[:] = data['PERR_'+name][:,1+8]
            newtab.data.field('Q_INNER')[:] = data[name+'FIT'][:,3+8]
            newtab.data.field('N_INNER')[:] = data[name+'FIT'][:,2+8]
            newtab.data.field('R_OUTER')[:] = data[name+'FIT'][:,1]
            newtab.data.field('R_OUTER_ERR')[:] = data['PERR_'+name][:,1]
            newtab.data.field('Q_OUTER')[:] = data[name+'FIT'][:,3]
            newtab.data.field('N_OUTER')[:] = data[name+'FIT'][:,2]
            newtab.data.field('INNER_FLUX_FRAC')[:] = data['FLUX_RATIO_'+name]
            newtab.data.field('R_TOT')[:] = data['REFF_'+name]
            newtab.data.field('R_TOT_ERR')[:] = np.sqrt(
                data['FLUX_RATIO_'+name]**2 *data['PERR_'+name][:,1+8]**2 + 
                (1-data['FLUX_RATIO_'+name])**2*data['PERR_'+name][:,1]**2 + 
                data['FLUX_RATIO_'+name]*(1-data['FLUX_RATIO_'+name])* 
                data['COVAR_'+name][:,1+(9*16)])
            newtab.data.field('TOTAL_FLUX')[:] \
                = (totalflux(data[name+'FIT'][:,0], 
                             data[name+'FIT'][:,1],
                             data[name+'FIT'][:,2],
                             data[name+'FIT'][:,3]) + 
                   totalflux(data[name+'FIT'][:,0+8], 
                             data[name+'FIT'][:,1+8],
                             data[name+'FIT'][:,2+8],
                             data[name+'FIT'][:,3+8]))
            newtab.data.field('I0_OUTER')[:] = data[name+'FIT'][:,0]*\
                np.exp(data[name+'FIT'][:,2]*np.exp(0.6950-0.1789/
                                                    data[name+'FIT'][:,2]))
            newtab.data.field('I0_INNER')[:] = data[name+'FIT'][:,0+8]*\
                np.exp(data[name+'FIT'][:,2+8]*np.exp(0.6950-0.1789/
                                                      data[name+'FIT'][:,2+8]))

        if not args.new:
            newtab.data.field('Q_TOT')[:] = data['SERSICFIT'][:,3]
        newtab.data.field('MAD_RATIO')[:] = (data['MAD_'+name+'_MASK'][:]/data['MAD_SKY'][:])
        newtab.data.field('RED_CHI2')[:] = data['CHISQ_'+name]
        newtab.data.field('STATUS')[:] = data['MPFIT_STATUS'][:,innums[n]]
        newtab.writeto(args.path+outfiles[n]+'.fits', clobber=True)

    for i in range(len(data.names)):
        print data.names[i], data.formats[i]

    return 0


if __name__=='__main__':
    main()
