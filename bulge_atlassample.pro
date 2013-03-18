;____________________________________
;
; Claire Lackner
;  Jan 14 2010
;
;  bulgeatlassample.pro
;
;  fits a bulge+disk to the atlas sample galaxies
;
; 
;____________________________________


;uses atlassample galaxies from atlasssampleDR7
PRO bulge_atlassample, filename, start, last, outputdir, imagedir, $
                       nosersic=nosersic, redo=redo, residuals=residuals, $
                       crop=crop
print, 'cpus: ',!CPU.TPOOL_NTHREADS     

gals = mrdfits(filename,1,range=[start,last-1])


;output entry for each object, the size of the array is for the bands,
;the default is five
if not keyword_set(redo) then begin
   output_entry = create_struct( gals[0], ['DSERSICFIT', $
                                           'SERSICFIT', $
                                           'DVCFIT', $
                                           'DDVCFIT', $
                                           'CHISQ_DSERSIC', $
                                           'CHISQ_SERSIC', $
                                           'CHISQ_DVC',$
                                           'CHISQ_DDVC', $
                                           'COVAR_DSERSIC', $
                                           'COVAR_SERSIC', $
                                           'COVAR_DVC', $
                                           'COVAR_DDVC', $
                                           'PERR_DSERSIC',$
                                           'PERR_SERSIC',$
                                           'PERR_DVC', $
                                           'PERR_DDVC', $
                                           'DOF_DSERSIC', $
                                           'DOF_SERSIC', $
                                           'DOF_DVC', $
                                           'DOF_DDVC', $
                                           'MPFIT_STATUS', $
                                           'XCROP', 'YCROP',$
                                           'XLEN', 'YLEN',$
                                           'MAD_SKY', $
                                           'MAD_DSERSIC', $
                                           'MAD_DSERSIC_MASK', $
                                           'MAD_SERSIC', $
                                           'MAD_SERSIC_MASK', $
                                           'MAD_DVC', $
                                           'MAD_DVC_MASK',$
                                           'MAD_DDVC', $
                                           'MAD_DDVC_MASK',$
                                           'SKY_DSERSIC', $
                                           'SKY_DSERSIC_ERR', $
                                           'SKY_DSERSIC_COVAR', $
                                           'SKY_SERSIC', $
                                           'SKY_SERSIC_ERR', $
                                           'SKY_SERSIC_COVAR', $
                                           'SKY_DVC', $
                                           'SKY_DVC_ERR', $
                                           'SKY_DVC_COVAR', $
                                           'SKY_DDVC', $
                                           'SKY_DDVC_ERR', $
                                           'SKY_DDVC_COVAR', $
                                           'FLUX_RATIO_DSERSIC', $
                                           'REFF_DSERSIC', $
                                           'FLUX_RATIO_DDVC', $
                                           'REFF_DDVC'], $
                                 dblarr(16), dblarr(8), $
                                 dblarr(8), dblarr(16), 0.0D, $
                                 0.0D, 0.0D, 0.0D, dblarr(16,16), $
                                 dblarr(8,8), dblarr(8,8), $
                                 dblarr(16,16), $
                                 dblarr(16), dblarr(8), $
                                 dblarr(8), dblarr(16), $
                                 0.0D, 0.0D, 0.0D, 0.0D, $
                                 dblarr(4), $
                                 0L, 0L, 0L, 0L, $
                                 0.0D, 0.0D, 0.0D, 0.0D, 0.0D, 0.0D, $
                                 0.0D, 0.0D, 0.0D, $
                                 0.0D, 0.0D, dblarr(17), $
                                 0.0D, 0.0D, dblarr(9), $
                                 0.0D, 0.0D, dblarr(9), $
                                 0.0D, 0.0D, dblarr(17), $
                                 0.0D, 0.0D, 0.0D, 0.0D )

endif else output_entry = create_struct(gals[0])

;insert known values into output strucutre array
output = replicate( output_entry, n_elements(gals) )
struct_assign, gals, output
;output.BANDS = (where(notignore ne 0))#(intarr(n_elements(gals))+1)

;file for writing out times
;openw, 1, 'timing.dat'

;loop over each image
for i=0L, n_elements(gals)-1L do begin

    print, gals[i].IDENT, ': ', i+1, ' of ', n_elements(gals)
    t1=systime(1)
        
    data = get_imivarpsf(gals[i].IDENT, strtrim(gals[i].FILENAME), imagedir)

    origsize = size(data.image, /dimensions )
    if n_elements(origsize) eq 1 then continue
    xcrop = 0
    ycrop = 0
    origdata=data
    if keyword_set(crop) then begin
       data = cropimage_2(data.image, data.ivar, data.psf, $
                          data.mask, $
                          x0=xcrop, y0=ycrop)
    endif else begin
       x0 = 0
       y0 = 0
    endelse
    
    imsize=size(data.image, /dimensions)
    if n_elements(imsize) lt 2 then begin
       print, '   cropped size: 0, skipping'
       continue
    endif

    print, '   cropped size:', imsize[0],'x',imsize[1]
    output[i].XLEN=imsize[0]
    output[i].YLEN=imsize[1]
    output[i].XCROP=xcrop
    output[i].YCROP=ycrop
        
    makegrid, imsize[0], imsize[1], x, y
    psf_fft = make_psf_fft(data.psf, data.image)

                                        ;fit dVc profile
    diskbulgefit,nod,params,data.image,psf_fft,data.ivar,$
      cs,covar,errors,fitstat,dof, sky, bulgeSersic=4.0, $
      nodisk=1, /free_sky, _EXTRA={Reff:0.2*sqrt(total(data.mask))}

    params[5] += xcrop
    params[6] += ycrop
    output[i].DVCFIT = params
    output[i].CHISQ_DVC = cs
    output[i].PERR_DVC[0:7] = errors[0:7]
    output[i].COVAR_DVC[0:7,0:7] = covar[0:7,0:7]          
    output[i].MPFIT_STATUS[2] = fitstat
    output[i].DOF_DVC = dof
    output[i].SKY_DVC = sky
    output[i].SKY_DVC_ERR = errors[8]
    output[i].SKY_DVC_COVAR = covar[8,0:8]

   
    if( not keyword_set(nosersic) ) then begin
        diskbulgefit, nod, params, data.image, psf_fft, data.ivar, $
          cs, covar, errors, fitstat, $
          dof, sky, bulgeSersic=4.0, $ ;4.0*fracDev + 1.0, $ 
          nodisk=1, freebulge=1, /free_sky, $
          _EXTRA={Reff:params[1], q:params[3], phi:params[7]}
                   
        params[5] += xcrop
        params[6] += ycrop
        output[i].SERSICFIT[*] = params
        output[i].CHISQ_SERSIC = cs
        output[i].PERR_SERSIC[0:7] = errors[0:7]
        output[i].COVAR_SERSIC[0:7,0:7] = covar[0:7,0:7]          
        output[i].MPFIT_STATUS[1] = fitstat
        output[i].DOF_SERSIC = dof
        output[i].SKY_SERSIC = sky
        output[i].SKY_SERSIC_ERR = errors[8]
        output[i].SKY_SERSIC_COVAR = covar[8,0:8]
     endif

    diskbulgefit, diskparam, bulgeparam, data.image, psf_fft, data.ivar, $
                  chsqds, covards, errds, stat, dofds, sky, bulgeSersic=4.0D, $
                  /free_sky, /freebulge, /freedisk, $
                  _EXTRA={Reff:params[1], q:params[3], phi:params[7], $
                          fracdev:0.8}
    output[i].MPFIT_STATUS[0] = stat
    output[i].SKY_DSERSIC = sky
    output[i].SKY_DSERSIC_ERR = errds[16]
    output[i].SKY_DSERSIC_COVAR = covards[16,0:16]
    diskparam[5:6] += [xcrop, ycrop]
    bulgeparam[5:6] += [xcrop, ycrop]
    output[i].PERR_DSERSIC = errds[0:15]
    output[i].COVAR_DSERSIC = covards[0:15,0:15]
    output[i].CHISQ_DSERSIC = chsqds
    output[i].DOF_DSERSIC = dofds
    if diskparam[2] ge bulgeparam[2] then $
       output[i].DSERSICFIT = [bulgeparam[*], diskparam[*]] $
    else $
       output[i].DSERSICFIT = [diskparam[*], bulgeparam[*]]

    output[i].FLUX_RATIO_DSERSIC = bulgetotot(output[i].DSERSICFIT, /cutoff)

    diskbulgefit, diskparam, bulgeparam, data.image, psf_fft, data.ivar, $
                  chsqds, covards, errds, stat, dofds, sky, bulgeSersic=4.0D, $
                  /free_sky, diskSersic=4.0D, $
                  _EXTRA={Reff:params[1]*4.0, q:params[3], phi:params[7], $
                          fracdev:0.1}
    output[i].MPFIT_STATUS[0] = stat
    output[i].SKY_DDVC = sky
    output[i].SKY_DDVC_ERR = errds[16]
    output[i].SKY_DDVC_COVAR = covards[16,0:16]
    diskparam[5:6] += [xcrop, ycrop]
    bulgeparam[5:6] += [xcrop, ycrop]
    output[i].PERR_DDVC = errds[0:15]
    output[i].COVAR_DDVC = covards[0:15,0:15]
    output[i].CHISQ_DDVC = chsqds
    output[i].DOF_DDVC = dofds
    if diskparam[1] lt bulgeparam[1] then $
       output[i].DDVCFIT = [bulgeparam[*], diskparam[*]] $
    else $
       output[i].DDVCFIT = [diskparam[*], bulgeparam[*]]

    output[i].FLUX_RATIO_DDVC = bulgetotot(output[i].DDVCFIT, /cutoff)

    print, 'timed: ',systime(1)-t1


    keep = where(origdata.mask,complement=notkeep)
    output[i].MAD_SKY = median(abs(origdata.image[notkeep]))


    ;put the model images in a file
    if keyword_set(residuals) then begin
        modName = string( output[i].IDENT, format='("models/M",i09,".fits")')
        modName = outputdir+modName
    endif
    makegrid, origsize[0], origsize[1], x, y
        
    model = pixelfluxpsf(x,y,$
                         [output[i].DSERSICFIT,output[i].SKY_DSERSIC], $
                         _EXTRA={cutoff:1,psfImage:data.psf})
        
    output[i].MAD_DSERSIC = median(abs(origdata.image-model))
    output[i].MAD_DSERSIC_MASK = median(abs((origdata.image-model)[keep]))
    half_flux = 0.5*totalsersicflux(output[i].DSERSICFIT, /cutoff)
    output[i].REFF_DSERSIC = $
       fluxfraction_radius(model-output[i].SKY_DSERSIC, $
                           fluxlevel=half_flux, $
                           q=output[i].SERSICFIT[3], $
                           phi=output[i].SERSICFIT[7], $
                           x0=output[i].DSERSICFIT[5], $
                           y0=output[i].DSERSICFIT[6], $
                           guess=[0.0, max([output[i].DSERSICFIT[9], $
                                            output[i].DSERSICFIT[1]])])

    if keyword_set(residuals) then begin
       mwrfits, model, modName, /create
    endif

    model = pixelfluxpsf(x,y,$
                         [output[i].SERSICFIT,output[i].SKY_SERSIC], $
                         _EXTRA={cutoff:1,psfImage:data.psf})
        
    output[i].MAD_SERSIC = median(abs(origdata.image-model))
    output[i].MAD_SERSIC_MASK = median(abs((origdata.image-model)[keep]))

    if keyword_set(residuals) then begin
       mwrfits, model, modName
    endif

    model = pixelfluxpsf(x,y,$
                         [output[i].DVCFIT,output[i].SKY_DVC], $
                         _EXTRA={cutoff:1,psfImage:data.psf})
        
    output[i].MAD_DVC = median(abs(origdata.image-model))
    output[i].MAD_DVC_MASK = median(abs((origdata.image-model)[keep]))

    if keyword_set(residuals) then begin
       mwrfits, model, modName
    endif

    model = pixelfluxpsf(x,y,$
                         [output[i].DDVCFIT,output[i].SKY_DDVC], $
                         _EXTRA={cutoff:1,psfImage:data.psf})
        
    output[i].MAD_DDVC = median(abs(origdata.image-model))
    output[i].MAD_DDVC_MASK = median(abs((origdata.image-model)[keep]))
    half_flux = 0.5*totalsersicflux(output[i].DDVCFIT, /cutoff)
    output[i].REFF_DDVC = $
       fluxfraction_radius(model-output[i].SKY_DDVC, $
                           fluxlevel=half_flux, $
                           q=output[i].DVCFIT[3], $
                           phi=output[i].DVCFIT[7], $
                           x0=output[i].DDVCFIT[5], $
                           y0=output[i].DDVCFIT[6], $
                           guess=[0.0, max([output[i].DDVCFIT[9], $
                                            output[i].DDVCFIT[1]])])
    if keyword_set(residuals) then begin
       mwrfits, model, modName
    endif


endfor
;close,1
;write output fits file: name of file=
append = string( start, last-1, format='("RAWFIT",i05,".",i05,".fits")')
outname = outputdir+append


mwrfits, output, outname, /create
print, 'writing output file: '+outname

END
