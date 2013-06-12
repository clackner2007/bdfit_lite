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
import matplotlib.cm
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
import matplotlib.figure as figure
from matplotlib.patches import Ellipse
import pyfits
from scipy.stats import scoreatpercentile
from ProfilePlots import getProfile

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

    for i in range(len(data)): #range(min(50,len(data))):
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
            modelImage = pyfits.open(args[2]+'M{0:09d}.fits'.\
                     format(int(data.IDENT[i])))[num[model]].data
            maskImage = pyfits.open(args[1]+'masks/{0}_mask.fits'.
                                    format(int(data.IDENT[i])))[0].data
        except IOError:
            s+= "<td>couldn't open image files for "+repr(data.IDENT[i])+"</td>\n</tr>\n"
            continue

        imX, imY = realImage.shape
        x0 = data[i]['XCROP']
        y0 = data[i]['YCROP']
        imX = data[i]['XLEN']
        imY = data[i]['YLEN']
        rs=3
        offset = min(data[i]['SERSICFIT'][1]*rs, 300)
        x0 = max(data[i]['SERSICFIT'][5] - offset, x0)
        y0 = max(data[i]['SERSICFIT'][6] - offset, y0)
        imX = min(imX, data[i]['SERSICFIT'][5] + offset - x0)
        imY = min(imY, data[i]['SERSICFIT'][6] + offset - y0)

        realImage[maskImage==1] = -1000
        showImage = copy.copy(realImage)
#    showImage[np.abs(ivarImage) < 1.e-13] = 0.0
        showImage=showImage[y0:y0+imY-1,x0:x0+imX-1]
        origImage = copy.copy(modelImage)
        modelImage[maskImage==1] = -1000
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
                   interpolation='none',
                   vmin=scoreatpercentile(np.arcsinh(showImage[showImage>-999].
                                                     flatten()), 5),
                   vmax=scoreatpercentile(np.arcsinh(showImage).flatten(),99.9))
        ax1.tick_params(labelsize='xx-small')

        ax2=fig.add_subplot(142, sharex=ax1, sharey=ax1)
        diff = (showImage-modelImage)
        diff[modelImage==0] = 0
        ax2.imshow(diff, aspect='equal',
                   interpolation='none',
                   vmin=scoreatpercentile(diff.flatten(), 0.5),
                   vmax=scoreatpercentile(diff.flatten(), 99.5))
        ax2.tick_params(labelleft='off',labelbottom='off')
        ax3=fig.add_subplot(143, sharex=ax1, sharey=ax1)
        plot_mini_model(ax3,
            bulgereff,diskreff,
            bulgephi,diskphi,
            bulgeq, diskq,
            bulgecenX-x0, bulgecenY-y0,
            diskcenX-x0, diskcenY-y0,
            imX, imY, btt)
        ax4=fig.add_subplot(144)
        realImage=np.asarray(realImage, dtype=np.float64)
        origImage=np.asarray(origImage, dtype=np.float64)
        realImage[maskImage==1] = -1000

        ax4.imshow(np.arcsinh(realImage),aspect='equal', 
                   interpolation='none',
                   vmin=scoreatpercentile(np.arcsinh(realImage
                                                     [realImage > 
                                                      -999]).flatten(), 0.1),
                   vmax=scoreatpercentile(np.arcsinh(showImage).flatten(),
                                          99.9))
        ax4.tick_params(labelsize='xx-small', labelleft='off', labelright='on')

        #improf = getProfile(realImage, data[i]['SERSICFIT'][1],
        #                    data[i]['SERSICFIT'][3],
        #                    data[i]['SERSICFIT'][7],
        #                    data[i]['SERSICFIT'][5],
        #                    data[i]['SERSICFIT'][6])
        #modprof = getProfile(origImage, data[i]['SERSICFIT'][1],
        #                     data[i]['SERSICFIT'][3],
        #                     data[i]['SERSICFIT'][7],
        #                     data[i]['SERSICFIT'][5],
        #                     data[i]['SERSICFIT'][6])
        #ax4.errorbar(improf.rad, improf.mnflux,
        #             yerr=improf.stdflux,
        #             fmt='o', mec='k', c='k')
        #ax4.plot(modprof.rad, modprof.mnflux, marker='None', 
        #         c='g', ls='solid')
        #ax4.tick_params(labelright='on', labelleft='off',
        #               labelbottom='on', labelsize='xx-small')
        #ax4.set_yscale('log')
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
