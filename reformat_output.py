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
        #wrong for everything but deVauc.
        flux *= 0.933666
    


def main():

    parser = argparse.ArgumentParser(description=__doc__,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("infile", help='input file')
    args = parser.parse_args()
    
    data = pyfits.open(args.infile)[1].data

    #make output files for each type of fit
    outcols = ['IDENT', 'R_INNER', 'N_INNER', 'R_OUTER', 'N_OUTER',
               'INNER_FLUX_FRAC', 'R_TOT', 'STATUS', 'Q_INNER',
               'Q_OUTER', 'Q_TOT', 'R_INNER_ERR', 'R_OUTER_ERR', 
               'R_TOT_ERR', 'TOTAL_FLUX']
    outtype = ['J', 'D', 'D', 'D', 'D', 'D', 'D', 'J', 'D', 'D', 'D',
               'D', 'D', 'D', 'D']
    

    outfiles = ['sersic', 'dblsersic', 'dbldvc', 'expsersic', 'dvc']
    innames = ['SERSIC', 'DSERSIC', 'DDVC', 'EXPSERSIC', 'DVC']
    
    nentry = len(data)
    
    for n, name in enumerate(innames):
        newtab = makeCols(outcols, outtype, nentry)
        print len(data['IDENT']), len(newtab.data['IDENT'])
        newtab.data.field('IDENT')[:] = data['IDENT']

        if (name=='SERSIC') or (name=='DVC'):
            newtab.data.field('R_INNER')[:] = data[name+'FIT'][:,1]
            newtab.data.field('R_INNER_ERR')[:] = data['PERR_'+name][:,1]
            newtab.data.field('Q_INNER')[:] = data[name+'FIT'][:,3]
            newtab.data.field('N_INNER')[:] = data[name+'FIT'][:,2]
            newtab.data.field('INNER_FLUX_FRAC')[:] = 1.0
            newtab.data.field('R_TOT')[:] = newtab.data.field('R_INNER')
            newtab.data.field('R_TOT_ERR')[:] = data['PERR_'+name][:,1]
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
            newtab.data.field('R_TOT_ERR')[:] = np.sqrt(\
                data['FLUX_RATIO_'+name]**2 *data['PERR_'+name][:,1+8]**2 + \
                (1-data['FLUX_RATIO_'+name])**2*data['PERR_'+name][:,1]**2 + \
                 data['FLUX_RATIO_'+name]*(1-data['FLUX_RATIO_'+name])* \
                 data['COVAR_'+name][:,1+(9*16)])
            newtab.data.field('TOTAL_FLUX')[:] \
                = totalflux(data[name+'FIT'][:,0], 
                            data[name+'FIT'][:,1],
                            data[name+'FIT'][:,2],
                            data[name+'FIT'][:,3]) + \
                            totalflux(data[name+'FIT'][:,0+8], 
                                      data[name+'FIT'][:,1+8],
                                      data[name+'FIT'][:,2+8],
                                      data[name+'FIT'][:,3+8])
                            

        newtab.data.field('Q_TOT')[:] = data['SERSICFIT'][:,3]
                
        newtab.writeto(outfiles[n]+'.fits', clobber=True)

    for i in range(len(data.names)):
        print data.names[i], data.formats[i]

    return 0


if __name__=='__main__':
    main()
