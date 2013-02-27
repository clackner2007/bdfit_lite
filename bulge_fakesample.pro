;____________________________________
;
; Claire Lackner
;  Jan 14 2010
;
; makes a fake galaxy image and fits it
; adjusts the leftover sky randomly, 50 times
;____________________________________



PRO bulge_fakesample, filename, start, last, nb=nb
;print, 'cpus: ',!CPU.TPOOL_NTHREADS

gals = mrdfits(filename,1,range=[start,last-1])

;number sky estimates
if not keyword_set(nb) then nb=1

folder = (strsplit(filename, '.fits', /extract, /regex ))[0]


output_entry = create_struct( gals[0], ['BULGEFIT', 'DISKFIT',$
                                        'SERSICFIT', 'CHISQ_BULGE', $
                                        'CHISQ_DISK', 'CHISQ_SERSIC', $
                                        'COVAR_BULGE', 'COVAR_DISK', $
                                        'COVAR_SERSIC', 'PERR_BULGE', $
                                        'PERR_DISK', 'PERR_SERSIC',$
                                        'MPFIT_STATUS', $
                                        'DOF_BULGE', 'DOF_DISK',$
                                        'DOF_SERSIC', 'DOF_DVC',$
                                        'DOF_EXP',$
                                        'EXPFIT', 'DVCFIT', $
                                        'CHISQ_EXP', 'CHISQ_DVC', $
                                        'PERR_EXP', 'PERR_DVC', $
                                        'COVAR_EXP', 'COVAR_DVC', $
                                        'BTT_BULGE', 'BTT_DISK', $
                                        'SKY_ERR', 'DISK_MAG', $
                                        'BULGE_MAG', 'MAG'], $
                              dblarr(14,nb), dblarr(14,nb), $
                              dblarr(8,nb), dblarr(nb), $
                              dblarr(nb), dblarr(nb), dblarr(14,14,nb), $
                              dblarr(14,14,nb), dblarr(8,8,nb), $
                              dblarr(14,nb), dblarr(14,nb), $
                              dblarr(8,nb), dblarr(5,nb),$
                              dblarr(nb),dblarr(nb),$
                              dblarr(nb),dblarr(nb),dblarr(nb),$
                              dblarr(8,nb), dblarr(8,nb),$
                              dblarr(nb), dblarr(nb), $
                              dblarr(8,nb),dblarr(8,nb),$
                              dblarr(8,8,nb), dblarr(8,8,nb), $
                              dblarr(nb), dblarr(nb), dblarr(nb), $
                              dblarr(nb), dblarr(nb), dblarr(nb) )
;insert known values into output strucutre array
output = replicate( output_entry, n_elements(gals) )
struct_assign, gals, output

seed=systime(1)

;sky leftover goes from -10,+10 DN
sky_err = randomu(seed, nb)*5 ;- 5
sky_err[0] = 0.

;loop over each image
for i=0L, n_elements(gals)-1L do begin

    print, i+1, ' of ', n_elements(gals)

    output[i].SKY_ERR=sky_err
    imsize=intarr(2) + $
      round(sqrt(gals[i].params[12]*gals[i].params[13]))*2
    print, imsize
    
    makegrid, imsize[0],imsize[1],x,y
    locpsf =  PSF_GAUSSIAN(npixel=51, $
                           fwhm=gals[i].psf_fwhm, /normalize)
    im_conv = pixelfluxpsf(x,y,gals[i].PARAMS, $
                           _EXTRA={psfImage:locpsf, cutoff:1})
    
    sky = gals[i].SKY_COUNTS
    ivar = 1./((im_conv+sky)/gals[i].GAIN)
    noise = randomn(seed,imsize[0],imsize[1])*sqrt(1./ivar)


    for sk=0, n_elements(sky_err) - 1 do begin

        ;add back in the error in the sky
        image = im_conv + noise + sky_err[sk]

        data = {image:image, ivar:ivar, psf:locpsf}
        xcrop = 0.0
        ycrop = 0.0
           

        reff = gals[i].PARAMS[7]*(randomn(seed)*0.1+1.0)
        rdisk = gals[i].PARAMS[1]*(randomn(seed)*0.1+1.0)
        reff = reff*gals[i].BULGE_TO_TOT + (1.-gals[i].BULGE_TO_TOT)*rdisk
        q = gals[i].PARAMS[9]*gals[i].BULGE_TO_TOT + $
          gals[i].PARAMS[3]*(1.-gals[i].BULGE_TO_TOT)
        phi=0.0
        fracdev=gals[i].BULGE_TO_TOT*(randomn(seed)*0.1+1.)

       ;  diskbulgefit, dbparams, bbparams, $
;           data.image, data.psf,$
;           data.ivar, chsqdb, covardb, errdb, stat, dofdb, $
;           bulgeSersic=4.0, $
;           _EXTRA={Reff:reff,$
;                   q:q,$
;                   phi:phi,$
;                   fracdev:fracdev}
        
;         output[i].MPFIT_STATUS[0,sk]=stat
        
        diskbulgefit, doparams, diparams, $
          data.image, data.psf, $
          data.ivar, chsqdd, covardd, errdd, stat, $
          dofdd, bulgeSersic=1.0, $
          _EXTRA={Reff:reff,$
                  q:q,$
                  phi:phi,$
                  fracdev:fracdev}
        
        output[i].MPFIT_STATUS[1,sk]=stat

        ;flip n=1 disk around, doesn't fix covariance
        if doparams[1] lt diparams[1] then begin
            temp = doparams
            doparams = diparams
            diparams = temp

            temp = errdd[0:5]
            errdd[0:5] = errdd[6:11]
            errdd[6:11] = temp

            temp = covardd[0:5,*]
            covardd[0:5,*] = covardd[6:11,*]
            covardd[6:11,*] = temp
            temp = covardd[*,0:5]
            covardd[*,0:5] = covardd[*,6:11]
            covardd[*,6:11] = temp
        endif
        
        diskdisk = [doparams[0:4],doparams[7],$
                    diparams[0:4],diparams[7],$
                    diparams[5:6]+[xcrop,ycrop]]
       ; dvcdisk = [dbparams[0:4],dbparams[7],$
       ;            bbparams[0:4],bbparams[7],$
       ;            bbparams[5:6]+[xcrop,ycrop]]
        
        output[i].DISKFIT[*,sk] = diskdisk
       ; output[i].BULGEFIT[*,sk] = dvcdisk
        ;output[i].PERR_BULGE[*,sk] = errdb
        output[i].PERR_DISK[*,sk] = errdd
        output[i].COVAR_DISK[*,*,sk] = covardd
        ;output[i].COVAR_BULGE[*,*,sk] = covardb
        ;output[i].CHISQ_BULGE[sk] = chsqdb
        output[i].CHISQ_DISK[sk] = chsqdd
        ;output[i].DOF_BULGE[sk] = dofdb
        output[i].DOF_DISK[sk] = dofdd
        ;output[i].BTT_BULGE[sk] = bulgetotot(dvcdisk)
        output[i].BTT_DISK[sk] = bulgetotot(diskdisk)

        flux = $
          totalsersicflux(diskdisk,cutoff=1)*output[i].NMGYPERCOUNT*$
          1.e-9
        output[i].MAG[sk] = sdss_asinh_mag(flux,filter=3)
        output[i].BULGE_MAG[sk] = sdss_asinh_mag(flux*output[i].BTT_DISK[sk],filter=3)
        output[i].DISK_MAG[sk] = $
          sdss_asinh_mag(flux*(1.-output[i].BTT_DISK[sk]),$
                         filter=3)
;         diskbulgefit, nod, params, data.image, data.psf, data.ivar, $
;           cs, covar, errors, fitstat, dof, bulgeSersic=4.0, $ 
;           nodisk=1, freebulge=1, $
;           _EXTRA={Reff:reff,$
;                   q:q,$
;                   phi:phi}

        
;         params[5] += xcrop
;         params[6] += ycrop
;         output[i].SERSICFIT[*,sk] = params
;         output[i].CHISQ_SERSIC[sk] = cs
;         output[i].PERR_SERSIC[*,sk] = errors
;         output[i].COVAR_SERSIC[*,*,sk] = covar           
;         output[i].MPFIT_STATUS[2,sk] = fitstat
;         output[i].DOF_SERSIC[sk] = dof
        
;         diskbulgefit,nod,params,data.image,data.psf,data.ivar,$
;           cs,covar,errors,fitstat,dof,bulgeSersic=4.0, $
;           nodisk=1, $
;           _EXTRA={Reff:reff,$
;                   q:q,$
;                   phi:phi}
        
;         params[5] += xcrop
;         params[6] += ycrop
;         output[i].DVCFIT[*,sk] = params
;         output[i].CHISQ_DVC[sk] = cs
;         output[i].PERR_DVC[*,sk] = errors
;         output[i].COVAR_DVC[*,*,sk] = covar           
;         output[i].MPFIT_STATUS[3,sk] = fitstat
;         output[i].DOF_DVC[sk] = dof
        
        diskbulgefit,nod,params,data.image,data.psf,data.ivar,$
          cs,covar,errors,fitstat,dof,bulgeSersic=1.0, $
          nodisk=1, $
          _EXTRA={Reff:rdisk,$
                  q:q,$
                  phi:phi}
        

        params[5] += xcrop
        params[6] += ycrop
        output[i].EXPFIT[*,sk] = params
        output[i].CHISQ_EXP[sk] = cs
        output[i].PERR_EXP[*,sk] = errors
        output[i].COVAR_EXP[*,*,sk] = covar           
        output[i].MPFIT_STATUS[4,sk] = fitstat
        output[i].DOF_EXP[sk] = dof
    endfor
endfor
    
;write output fits file: name of file
;sp=strsplit(folder,'/',/extract)
;outfold=sp[n_elements(sp)-1]
outfold=folder


append = string( start, last-1, format='("RAWFIT",i05,".",i05,".fits")') 

outname = outfold + '/fakePOS_' + append

mwrfits, output, outname, /silent, /create
print, 'writing output file: ',outname

END
