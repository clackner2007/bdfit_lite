#!/usr/bin/env python

##############################
# make galaxy images for web display
#
##############################
"""
%prog fitter_raw_output image_folder model_folder output_image/folder
makes the web output for the cosmos images for a given model
"""

import sys
import os
import re, copy
import optparse
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
import matplotlib.figure as figure
from matplotlib.patches import Ellipse
import pyfits
from scipy.stats import scoreatpercentile

##plot stylized model
def plot_mini_model(ax,
            bulgereff,diskreff,
            bulgephi,diskphi,
            bulgeq, diskq,
            bulgexcen, bulgeycen,
            diskxcen, diskycen,
            x_len, y_len, btt):

    if diskreff > 0:
        ax.add_artist((Ellipse((diskxcen,diskycen),diskreff,
                               diskq*diskreff,
                               diskphi*180./np.pi,ec='b',fc='b',
                               linewidth=0.0,alpha=0.4*(1.0-btt)+0.1)))
        ax.add_artist((Ellipse((diskxcen,diskycen),diskreff,
                              diskq*diskreff,
                              diskphi*180./np.pi,ec='b',fc='None',
                              linewidth=0.5,alpha=1.0)))

    if bulgereff > 0:
        ax.add_artist((Ellipse((bulgexcen,bulgeycen),bulgereff,
                               bulgeq*bulgereff,
                               bulgephi*180./np.pi,ec='r',fc='r',
                               linewidth=0.0,alpha=(btt)*0.4+0.1)))
        ax.add_artist((Ellipse((bulgexcen,bulgeycen),bulgereff,
                   bulgeq*bulgereff,
                   bulgephi*180./np.pi,ec='r',fc='None',
                   linewidth=0.5)))

    ax.set_aspect('equal')
    ax.set_xlim((0,x_len))
    ax.set_ylim((y_len,0))
    ax.tick_params(labelleft='off',labelright='on',labelsize='xx-small')


##



def main():
    parser = optparse.OptionParser(usage=__doc__)
    parser.add_option('-m', '--model',
                      help='model to use, bulge, disk, exp, dvc, or sersic',
                      dest='model', default='sersic')

    opts, args = parser.parse_args()

    np.seterr(all='ignore')
    if len(args) != 4:
        parser.print_help()
        sys.exit(1)


    #read the data
    data=pyfits.open(args[0])[1].data
    num={'DSERSIC':0,'DDVC':3,'EXPSERSIC':4,'DVC':2,'SERSIC':1}

    model = str.upper(opts.model)
    if model not in num.keys():
        model='SERSIC'

    #write the n_b=4 b+d page
    s  = "<html>\n"
    s += "<body>\n"
    s += "<table border=\"1\">\n"

    for i in range(min(50,len(data))):
        s += "<tr>\n"
        s2 = "<table>\n"
        keys=['IDENT',#'ALPHA_J2000','DELTA_J2000',
            'CHISQ_%s'%model,
          #'FRACDEV', #'ZEST_TYPE',
            'MAD_SKY']
        if model in ['DSERSIC', 'DDVC', 'EXPSERSIC']:
            keys.append('REFF_'+model)
            keys.append('FLUX_RATIO_'+model)

        keys.append('MAD_%s_MASK'%model)
    #keys=('IDENT','ALPHA_J2000','DELTA_J2000','CHISQ_BULGE',
    #      'FRACDEV')
        for k in keys:
            v = data[i][k]
            s2 +=  "<tr><td>%s:</td><td>%s</td></tr>\n" % (k, str(v))

        s2 +=  "<tr><td>%s</td><td>%s</td></tr>\n" % ('SERSIC N:', str(data[i]['SERSICFIT'][2]))
        if model=='DSERSIC':
            s2 +=  "<tr><td>%s</td><td>%.2f,%.2f</td></tr>\n" % ('DSERSIC NS:', data[i]['DSERSICFIT'][2],
                                                                 data[i]['DSERSICFIT'][10])
        if model=='EXPSERSIC':
            s2 += "<tr><td>%s</td><td>%s</td></tr>\n" % ('EXPSERSIC N:', str(data[i]['EXPSERSICFIT'][10]))
#    s2 +=  "<tr><td>%s</td><td>%s</td></tr>\n" % ('B/T:', str(data[i]['SERSICFIT'][2]))
        s2 += "</table>\n"
        s += "<td>%s</td>\n" % (s2)

    #get the images
        try:
            realImage = pyfits.open(args[1]+'images/{0}.fits'. \
                    format(data.FILENAME[i]))[0].data
            #ivarImage = pyfits.open(args[1]+'ivar/{0}.0_{1:f}_{2:f}.wht.mask.fits'.\
            #        format(data.IDENT[i],data.ALPHA_J2000[i],data.DELTA_J2000[i]))[0].data
            #psfImage = pyfits.open(args[1]+'psf/{0}.0_{1:f}_{2:f}.psf.fits.gz'.\
            #       format(data.IDENT[i],data.ALPHA_J2000[i],data.DELTA_J2000[i]))[0].data
            modelImage = pyfits.open(args[2]+'M{0:09d}.fits'.\
                     format(int(data.IDENT[i])))[num[model]].data
        except IOError:
            s+= "<td>couldn't open image files for "+repr(data.IDENT[i])+"</td>\n</tr>\n"
            continue

        imX, imY = realImage.shape
        x0 = data[i]['XCROP']
        y0 = data[i]['YCROP']
        imX = data[i]['XLEN']
        imY = data[i]['YLEN']
#    xx, yy = np.meshgrid(range(imX), range(imY))
#    x0 = np.min(xx[ivarImage > 1.e-13])
#    y0 = np.min(yy[ivarImage > 1.e-13])
#    imX = np.max(xx[ivarImage > 1.e-13]) - x0
#    imY = np.max(yy[ivarImage > 1.e-13]) - y0
#    del xx, yy

        showImage = realImage #copy.copy(realImage)
#    showImage[np.abs(ivarImage) < 1.e-13] = 0.0
        showImage=showImage[y0:y0+imY-1,x0:x0+imX-1]
        modelImage=modelImage[y0:y0+imY-1,x0:x0+imX-1]

        modelfit='%sFIT'%model
        if model in ('DDVC', 'EXPSERSIC', 'DSERSIC'):
            bulgereff=data[i][modelfit][9]
            bulgeq=data[i][modelfit][11]
            bulgephi=data[i][modelfit][15]
            diskreff=data[i][modelfit][1]
            diskq=data[i][modelfit][3]
            diskphi=data[i][modelfit][7]
            bulgecenX=data[i][modelfit][13]
            bulgecenY=data[i][modelfit][14]
            diskcenX=data[i][modelfit][5]
            diskcenY=data[i][modelfit][6]
            btt={'DDVC':data[i]['FLUX_RATIO_DDVC'],
                 'EXPSERSIC':data[i]['FLUX_RATIO_EXPSERSIC'],
                'DSERSIC':data[i]['FLUX_RATIO_DSERSIC']}[model]
        else:
            bulgereff=data[i][modelfit][1]
            bulgeq=data[i][modelfit][3]
            bulgephi=data[i][modelfit][7]
            diskreff=0.0
            diskq=0.0
            diskphi=0.0
            bulgecenX=data[i][modelfit][5]
            bulgecenY=data[i][modelfit][6]
            diskcenX=0.0
            diskcenY=0.0
            btt=1.0
      #bulgereff=data[i]['SERSICFIT'][1]
      #bulgeq=data[i]['SERSICFIT'][3]
      #bulgephi=data[i]['SERSICFIT'][7]
      #bulgecenX=data[i]['SERSICFIT'][5]
      #bulgecenY=data[i]['SERSICFIT'][6]
      #diskcenX=0.0
      #diskcenY=0.0
      #diskreff=0.0
      #diskq=0.0
      #diskphi=0.0
      #btt=1.0


        fig=figure.Figure((8.16,imY*1.0/imX*2.0))
        fig.subplots_adjust(wspace=0,left=0.1,right=0.9,
                            top=1.0,bottom=0.0,hspace=0.0)
        canv = FigCanvas(fig)

        ax1=fig.add_subplot(141,frameon=True)
        ax1.imshow((np.arcsinh(showImage)),aspect='equal',
                   vmin=np.amin(np.arcsinh(showImage)),
                   vmax=scoreatpercentile(np.arcsinh(showImage).flatten(),99.9))
        ax1.tick_params(labelsize='xx-small')

        ax2=fig.add_subplot(142)
        ax2.imshow((showImage-modelImage),aspect='equal')
        ax2.tick_params(labelleft='off',labelbottom='off')
        ax3=fig.add_subplot(143)
        plot_mini_model(ax3,
            bulgereff,diskreff,
            bulgephi,diskphi,
            bulgeq, diskq,
            bulgecenX-x0, bulgecenY-y0,
            diskcenX-x0, diskcenY-y0,
            imX, imY, btt)
        #ax4=fig.add_subplot(144)
        #ax4.imshow(np.arcsinh(psfImage),aspect='equal')
        #ax4.tick_params(labelleft='off',labelbottom='off',labelright='on',
        #    labelsize='xx-small')
        #ax4.set_xlim(0,imX)
        #ax4.set_ylim(imY,0)

        fig.savefig('{2}/{1}fit/{0:09d}.png'.format(int(data.IDENT[i]),
             str.lower(model), args[3]))
        s += "<td><img src=\"{1}fit/{0:09d}.png\"></td>\n".format(int(data.IDENT[i]),
                                                            str.lower(model))
        s += "</tr>\n"


    s += "</table>\n"
    s += "</body>\n"
    s += "</html>\n"

    print s


if __name__=='__main__':
    main()
