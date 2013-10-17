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
PRO bulge_atlassample2, filename, start, last, outputdir, imagedir, $
                       nosersic=nosersic, redo=redo, residuals=residuals, $
                       crop=crop
print, 'cpus: ',!CPU.TPOOL_NTHREADS     

gals = mrdfits(filename,1,range=[start,last-1])
inputs = mrdfits('bcg_data/outputs/final_match.fits', 1, range=[start,last-1])

;output entry for each object, the size of the array is for the bands,
;the default is five
if not keyword_set(redo) then begin
   output_entry = create_struct( gals[0], ['DVCSERSICFIT', $
                                           'DVCEXPFIT', $
                                           'DDVCFIT', $
                                           'CHISQ_DVCSERSIC', $
                                           'CHISQ_DVCEXP',$
                                           'CHISQ_DDVC', $
                                           'COVAR_DVCSERSIC', $
                                           'COVAR_DVCEXP', $
                                           'COVAR_DDVC', $
                                           'PERR_DVCSERSIC',$
                                           'PERR_DVCEXP', $
                                           'PERR_DDVC', $
                                           'DOF_DVCSERSIC', $
                                           'DOF_DVCEXP', $
                                           'DOF_DDVC', $
                                           'MPFIT_STATUS', $
                                           'XCROP', 'YCROP',$
                                           'XLEN', 'YLEN',$
                                           'MAD_SKY', $
                                           'MAD_DVCSERSIC', $
                                           'MAD_DVCSERSIC_MASK', $
                                           'MAD_DVCEXP', $
                                           'MAD_DVCEXP_MASK',$
                                           'MAD_DDVC', $
                                           'MAD_DDVC_MASK',$
                                           'SKY_DVCSERSIC', $
                                           'SKY_DVCSERSIC_ERR', $
                                           'SKY_DVCSERSIC_COVAR', $
                                           'SKY_DVCEXP', $
                                           'SKY_DVCEXP_ERR', $
                                           'SKY_DVCEXP_COVAR', $
                                           'SKY_DDVC', $
                                           'SKY_DDVC_ERR', $
                                           'SKY_DDVC_COVAR', $
                                           'FLUX_RATIO_DVCSERSIC', $
                                           'REFF_DVCSERSIC', $
                                           'FLUX_RATIO_DDVC', $
                                           'REFF_DDVC', $
                                           'FLUX_RATIO_DVCEXP', $
                                           'REFF_DVCEXP'], $
                                 dblarr(16), dblarr(16), dblarr(16), $
                                 0.0D, 0.0D, 0.0D, $
                                 dblarr(16,16), $
                                 dblarr(16,16), dblarr(16,16), $
                                 dblarr(16), dblarr(16), dblarr(16), $
                                 0.0D, 0.0D, 0.0D, $
                                 dblarr(3), $
                                 0L, 0L, 0L, 0L, $
                                 0.0D, $
                                 0.0D, 0.0D, 0.0D, $
                                 0.0D, 0.0D, 0.0D, $
                                 0.0D, 0.0D, dblarr(17), $
                                 0.0D, 0.0D, dblarr(17), $
                                 0.0D, 0.0D, dblarr(17), $
                                 0.0D, 0.0D, 0.0D, 0.0D, 0.0D, 0.0D )

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

    ;get old fits info for initial conditions
    params=inputs[i].SERSICFIT
    reff_in = min([params[1]*1.2, sqrt(imsize[0]*imsize[1]*1.0)*0.8])
    
    diskbulgefit, diskparam, bulgeparam, data.image, psf_fft, data.ivar, $
                  chsqds, covards, errds, stat, dofds, sky, bulgeSersic=4.0D, $
                  /free_sky, diskSersic=1.0D, /freedisk, /rescale, $
                  _EXTRA={Reff:reff_in, q:params[3], phi:params[7], $
                          fracdev:0.5}
    output[i].MPFIT_STATUS[0] = stat
    output[i].SKY_DVCSERSIC = sky
    output[i].SKY_DVCSERSIC_ERR = errds[16]
    output[i].SKY_DVCSERSIC_COVAR = covards[16,0:16]
    diskparam[5:6] += [xcrop, ycrop]
    bulgeparam[5:6] += [xcrop, ycrop]
    output[i].PERR_DVCSERSIC = errds[0:15]
    output[i].COVAR_DVCSERSIC = covards[0:15,0:15]
    output[i].CHISQ_DVCSERSIC = chsqds
    output[i].DOF_DVCSERSIC = dofds
    if diskparam[2] ge bulgeparam[2] then $
       output[i].DVCSERSICFIT = [bulgeparam[*], diskparam[*]] $
    else $
       output[i].DVCSERSICFIT = [diskparam[*], bulgeparam[*]]

    output[i].FLUX_RATIO_DVCSERSIC = bulgetotot(output[i].DVCSERSICFIT, /cutoff)

    diskbulgefit, diskparam, bulgeparam, data.image, psf_fft, data.ivar, $
                  chsqds, covards, errds, stat, dofds, sky, bulgeSersic=4.0D, $
                  /free_sky, diskSersic=4.0D, $
                  _EXTRA={Reff:reff_in, $
                          q:params[3], phi:params[7], $
                          fracdev:0.5}
    output[i].MPFIT_STATUS[1] = stat
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

    diskbulgefit, diskparam, bulgeparam, data.image, psf_fft, data.ivar, $
                  chsqds, covards, errds, stat, dofds, sky, bulgeSersic=4.0D, $
                  /free_sky, diskSersic=1.0D, /rescale, $
                  _EXTRA={Reff:reff_in, $
                          q:params[3], phi:params[7], $
                          fracdev:0.5}
    output[i].MPFIT_STATUS[2] = stat
    output[i].SKY_DVCEXP = sky
    output[i].SKY_DVCEXP_ERR = errds[16]
    output[i].SKY_DVCEXP_COVAR = covards[16,0:16]
    diskparam[5:6] += [xcrop, ycrop]
    bulgeparam[5:6] += [xcrop, ycrop]
    output[i].PERR_DVCEXP = errds[0:15]
    output[i].COVAR_DVCEXP = covards[0:15,0:15]
    output[i].CHISQ_DVCEXP = chsqds
    output[i].DOF_DVCEXP = dofds
    
    output[i].DVCEXPFIT = [diskparam[*], bulgeparam[*]]

    output[i].FLUX_RATIO_DVCEXP = bulgetotot(output[i].DVCEXPFIT, /cutoff)

    print, 'timed: ',systime(1)-t1
    ;mwrfits, output, 'temp.fits', /create

    keep = where(origdata.mask,complement=notkeep)
    output[i].MAD_SKY = median(abs(origdata.image[notkeep]))


    ;put the model images in a file
    if keyword_set(residuals) then begin
        modName = string( output[i].IDENT, format='("models/M",i09,".fits")')
        modName = outputdir+modName
    endif
    makegrid, origsize[0], origsize[1], x, y
        
    model = pixelfluxpsf(x,y,$
                         [output[i].DVCSERSICFIT,output[i].SKY_DVCSERSIC], $
                         _EXTRA={cutoff:1,psfImage:data.psf})
        
    output[i].MAD_DVCSERSIC = median(abs(origdata.image-model))
    output[i].MAD_DVCSERSIC_MASK = median(abs((origdata.image-model)[keep]))
    half_flux = 0.5*totalsersicflux(output[i].DVCSERSICFIT, /cutoff)
    if output[i].MPFIT_STATUS[0] eq 0 then $
       output[i].REFF_DVCSERSIC = -1.0 else $
          output[i].REFF_DVCSERSIC = $
       fluxfraction_radius(model-output[i].SKY_DVCSERSIC, $
                           fluxlevel=half_flux, $
                           q=inputs[i].SERSICFIT[3], $
                           phi=inputs[i].SERSICFIT[7], $
                           x0=output[i].DVCSERSICFIT[5], $
                           y0=output[i].DVCSERSICFIT[6], $
                           guess=[0.0, max([output[i].DVCSERSICFIT[9], $
                                            output[i].DVCSERSICFIT[1]])])

    if keyword_set(residuals) then begin
       mwrfits, model, modName, /create
    endif

  
    model = pixelfluxpsf(x,y,$
                         [output[i].DDVCFIT,output[i].SKY_DDVC], $
                         _EXTRA={cutoff:1,psfImage:data.psf})
        
    output[i].MAD_DDVC = median(abs(origdata.image-model))
    output[i].MAD_DDVC_MASK = median(abs((origdata.image-model)[keep]))
    half_flux = 0.5*totalsersicflux(output[i].DDVCFIT, /cutoff)
    if output[i].MPFIT_STATUS[1] eq 0 then $
       output[i].REFF_DDVC = -1.0 else $
          output[i].REFF_DDVC = $
       fluxfraction_radius(model-output[i].SKY_DDVC, $
                           fluxlevel=half_flux, $
                           q=inputs[i].SERSICFIT[3], $
                           phi=inputs[i].SERSICFIT[7], $
                           x0=output[i].DDVCFIT[5], $
                           y0=output[i].DDVCFIT[6], $
                           guess=[0.0, max([output[i].DDVCFIT[9], $
                                            output[i].DDVCFIT[1]])])
    if keyword_set(residuals) then begin
       mwrfits, model, modName
    endif
    
    model = pixelfluxpsf(x,y,$
                         [output[i].DVCEXPFIT,output[i].SKY_DVCEXP], $
                         _EXTRA={cutoff:1,psfImage:data.psf})
        
    output[i].MAD_DVCEXP = median(abs(origdata.image-model))
    output[i].MAD_DVCEXP_MASK = median(abs((origdata.image-model)[keep]))
    half_flux = 0.5*totalsersicflux(output[i].DVCEXPFIT, /cutoff)
    if output[i].MPFIT_STATUS[2] eq 0 then $
       output[i].REFF_DVCEXP = -1.0 else $
          output[i].REFF_DVCEXP = $
       fluxfraction_radius(model-output[i].SKY_DVCEXP, $
                           fluxlevel=half_flux, $
                           q=inputs[i].SERSICFIT[3], $
                           phi=inputs[i].SERSICFIT[7], $
                           x0=output[i].DVCEXPFIT[5], $
                           y0=output[i].DVCEXPFIT[6], $
                           guess=[0.0, max([output[i].DVCEXPFIT[9], $
                                            output[i].DVCEXPFIT[1]])])
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
