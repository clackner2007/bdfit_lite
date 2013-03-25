#!usr/bin/env python
############################################
# Claire Lackner
# June 16 2010
#
# ProfilePlots.py
# given the atlas and model image(s), 
# this plots a 'radial profile' by taking bins of the model image 
#
#
############################################
"""
%prog
"""
import sys, os, re
import optparse
import math as m
#use cython version
from repixelate import repixel
import numpy as np


import pyfits
import matplotlib.pyplot as plt

class Ellipse:
    def __init__(self, r, q, x0, y0, phi):
	self.r = r
	self.x0 = x0
	self.y0 = y0
	self.phi = phi
	self.q = q
   
    def recenter(self,x,y):
	return ((x-self.x0)*np.cos(self.phi) + (y-self.y0)*np.sin(self.phi),
		(y-self.y0)*np.cos(self.phi) - (x-self.x0)*np.sin(self.phi))

    #retrun true if x,y point is in ellipse
    def is_inside( self, x, y ):
	(xt,yt)=self.recenter(x,y)
	inside=(xt<self.r)&(yt<self.r)
	small=np.where(inside)

	rads=np.sqrt(xt[small]*xt[small] + yt[small]*yt[small]/(self.q*self.q))
        inside[small] = rads < self.r
	return inside


##################################
# determine the total flux within an annulus
# also return the mean and the variance
def FluxAnnulus( image, r_in, r_out, q, phi, x0, y0 ):
    #repixel by 1
    rp= 1.
    image_64 = repixel(np.asarray(image), 1/rp )
    norm = np.sum(image)/np.sum(image_64)

    #image_64 = image
    x,y = np.meshgrid( np.arange(np.shape(image_64)[1]*1.0),
		       np.arange(np.shape(image_64)[0]*1.0) )
    ell1 = ~(Ellipse(r_in*rp, q, rp*x0, rp*y0, phi ).is_inside(x,y))
    ell2 = Ellipse(r_out*rp, q, rp*x0, rp*y0, phi).is_inside(x,y)
    annulus = ell1 & ell2
    if ~annulus.any():
	return (0.0,0.0,0.0,-1.0)
    totalFlux,meanFlux,varFlux,area=(0.0,0.0,0.,0.)
    totalFlux = np.sum(image_64[annulus])*norm
    area = (image_64[annulus].size/(rp*rp))
    meanFlux = np.mean(image_64[annulus]*norm)
    varFlux = np.var(image_64[annulus]*norm)
    if(totalFlux == 0.0): area = -1.0
    return (totalFlux, meanFlux, varFlux, area)
		       

def getProfile(image, reff, q, phi, x0, y0):

    step = min(max(0.05*reff,0.8),1.5)
    r=0.0
    numstep=18
    high = 4.*reff
    mult_factor = 10.**(np.log10(high-step)/numstep)
 
    profile = np.recarray((numstep,),dtype=[('rad',float),
					    ('mnflux',float),
					    ('stdflux',float),
					    ('sb', float)])
    for i in range(numstep):
	profile[i].rad = 0.5*(r+step)
	(tf, mf, vf, area) = FluxAnnulus(image, r, r+step, q, phi, x0, y0)
	profile[i].mnflux = mf
	profile[i].stdflux = m.sqrt(vf)
	profile[i].sb = tf/area
	r = r + step
	step *= mult_factor

    return profile
	
#####################
# for testing
def main():
    
    imagefile=pyfits.open('/u/clackner/work/bulges/output_06132010/S0000-00094/images/SDSSJ022601.26+010105.0-n4.fits')
    im2=pyfits.open('/u/clackner/work/bulges/output_06132010/S0000-00094/images/SDSSJ022601.26+010105.0-n1.fits')
    image=np.asarray(imagefile[0].data-imagefile[1].data-1000.0)
    model=np.asarray(imagefile[0].data-1000.0)
    bmod=np.asarray(imagefile[3].data,dtype=np.float64)
    dmod=np.asarray(model-bmod)
    mod2=im2[0].data -1000.0

    data=pyfits.open('/u/clackner/work/bulges/output_06132010/S0000-00094/totalRAW.fits')[1].data
    want=np.where(data.GALID == 'SDSSJ022601.26+010105.0')

    bulge = data[want].field('BULGEFIT')[:,7]
    disk = data[want].field('BULGEFIT')[:,1]
    qb = data[want].field('BULGEFIT')[:,9]
    qd = data[want].field('BULGEFIT')[:,3]
    phib=data[want].field('BULGEFIT')[:,11]
    phid=data[want].field('BULGEFIT')[:,5]
    xc=data[want].field('BULGEFIT')[:,12]
    yc=data[want].field('BULGEFIT')[:,13]

    if(bulge>disk):
	reff=bulge
	q=qb
	phi=phib
    else:
	reff=disk
	q=qd
	phi=phid
#    print reff, q, phi
    
    improf = getProfile(image, reff, q, phi, xc, yc )
    modprof = getProfile(model, reff,q,phi,xc,yc)
    mod2prof = getProfile(mod2, reff,q,phi,xc,yc)
    bmodprof = getProfile(bmod, reff,q,phi,xc,yc)
    dmodprof = getProfile(dmod, reff,q,phi,xc,yc)
 
    plt.xscale('log')
    plt.errorbar(modprof.rad,modprof.sb,yerr=None,marker='^',
		 mec='g',mfc='None',mew=1.4,label='model, n=4',c='g', ms=8.)
    plt.plot(bmodprof.rad, bmodprof.sb, marker='^', c='r', mec='r',label='bulge')
    plt.plot(dmodprof.rad, dmodprof.sb, marker='^', c='b', mec='b',label='disk')
    plt.errorbar(improf.rad,improf.sb,yerr=improf.stdflux, marker='.',
		 mfc='k',label='image',c='k')
    plt.ylim((0.0,np.max([improf.sb,modprof.sb])*1.02))
    plt.xlim((improf.rad[0]*0.9, improf.rad[-1]))
    plt.legend()
    plt.xlabel('radius in pixels')
    plt.ylabel('surface brightness in counts/pix^2')
    plt.show()

if __name__=='__main__':
    main()
