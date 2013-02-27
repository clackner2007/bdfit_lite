;____________________________________
;
; Claire Lackner
;  Jan 14 2010
;
;  change the s/n or the resolution of a  galaxy as it's being fit
;
; 
;____________________________________


;uses atlassample galaxies from atlasssampleDR7
;do the fits with changing sn/resolve for fake galaxies
PRO bulge_sn_resolve_fake, filename, start, last, outputdir, $
                           nosersic=nosersic, noexp=noexp, nodvc=nodvc, $
                           redo=redo, $
                           add_fwhm=add_fwhm, add_noise=add_noise
;pref_set, 'idl_cpu_tpool_nthreads', 0, /commit
print, 'cpus: ',!CPU.TPOOL_NTHREADS     

gals = mrdfits(filename,1,range=[start,last-1])

nb=1
notignore=[0,0,1,0,0]
;output entry for each object, the size of the array is for the bands,
;the default is five
if not keyword_set(redo) then begin
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
                                            'BANDS', 'ADD_NOISE', 'ADD_FWHM'], $
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
                                  intarr(nb),dblarr(nb),dblarr(nb) )
    
endif else output_entry = create_struct(gals[0], ['ADD_NOISE', 'ADD_FWHM'],$
                                        dblarr(nb),dblarr(nb))

;triple sky counts
gals.SKY_COUNTS *= 40.

;insert known values into output strucutre array
output = replicate( output_entry, n_elements(gals) )
struct_assign, gals, output
;output.BANDS = reform((where(notignore ne 0))#(intarr(n_elements(gals))+1),n_elements(gals),1)

;file for writing out times
;openw, 1, 'timing.dat'

;loop over each image
for i=0L, n_elements(gals)-1L do begin

    print, i+1, ' of ', n_elements(gals)
    gals[i].params[3] = 1.0
    output[i].params[3] = 1.0
    ;double the image size to make sure weird cutoff aren't happening.
    gals[i].params[5] *= 2.0
    gals[i].params[6] *= 2.0
    output[i].params[5] *= 2.0
    output[i].params[6] *= 2.0
    ;half the scale elgnth
    gals[i].params[1] *= 0.5
    output[i].params[1] *= 0.5

        
    ;get image from parameter file
    imsize=intarr(2) + round(sqrt(gals[i].params[5]*gals[i].params[6]))*2
    
    makegrid, imsize[0],imsize[1],x,y
    locpsf =  0.9*PSF_GAUSSIAN(npixel=51, $
                               fwhm=gals[i].psf_fwhm, /normalize) + $
      0.1*PSF_GAUSSIAN(npixel=51,fwhm=gals[i].psf_fwhm*2., /normalize)
    if keyword_set(add_fwhm) then begin
        locpsf = 0.9*PSF_GAUSSIAN(npixel=51, $
                                  fwhm=gals[i].psf_fwhm*add_fwhm, /normalize) + $
          0.1*PSF_GAUSSIAN(npixel=51,fwhm=gals[i].psf_fwhm*2.*add_fwhm, /normalize)
        output[i].ADD_FWHM[0] = add_fwhm
    endif
    im_conv = pixelfluxpsf(x,y,gals[i].PARAMS, $
                           _EXTRA={psfImage:locpsf, cutoff:0})
    sky = gals[i].SKY_COUNTS
    if keyword_set(add_noise) then begin
        sky *= add_noise^2
        output[i].ADD_NOISE = add_noise
    endif
    ivar = 1./((im_conv+sky)/gals[i].GAIN)
    noise = randomn(seed,imsize[0],imsize[1])*sqrt(1./ivar)
    image = im_conv + noise
    ivar = 1./((1.*image + sky)/gals[i].GAIN)
    gain = gals[i].GAIN

    
;                                         ;convolve the image with a larger psf
;     if keyword_set(add_fwhm) then begin
;         resize=1
;         delta_psffwhm = (add_fwhm-1.)*gals[i].PSF_FWHM

      
;         noise= ivar
;         noise[where(noise ne 0)] = 1./sqrt(noise[where(noise ne 0)])
        
;         psf_size=101
;         add_psf_correct = PSF_GAUSSIAN(npixel=psf_size, $
;                                        fwhm=delta_psffwhm, /normalize)
;         im_large = convolve(image, add_psf_correct $
;                                 )
;         psf_large = convolve(locpsf, add_psf_correct $
;                             )
;         noise_large = convolve(noise, add_psf_correct $
;                                    )
;         add_var = (1.-2.*alog(2.)/(!pi*(delta_psffwhm)^2))* $
;           sqrt((image+sky)/gain)
;         im_large += (randomn(seed, imsize[0], $
;                              imsize[1]))*sqrt(add_var)

  
;         image = im_large 
;         ;image[where(ivar eq 0)] = 0.
;         locpsf = psf_large         
;         noise2_large = noise_large^2 + add_var

                                        
;         ivar[where((noise2_large ne 0)and(ivar ne 0))] = $
;           1./noise2_large[where((noise2_large ne 0)and(ivar ne 0))]
;         output[i].ADD_FWHM[0] = add_fwhm
;     endif

;         ;increase the total amount noise by upping the sky by add_noise^2
;      if keyword_set(add_noise) then begin
;          sky_noise = sqrt(sky/gain)
        
;         seed = systime(1)
;         image[where(ivar ne 0)] += ((randomn(seed, imsize[0], $
;                                                  imsize[1]))[where(ivar ne 0)])*$
;           sky_noise * $           
;           sqrt(add_noise*add_noise-1.)
;         var = ivar
;         var[where(ivar ne 0)] = 1./ivar[where(ivar ne 0)] + $
;           sky_noise^2*(add_noise^2-1.)
;         ivar[where(ivar ne 0)] = 1./var[where(ivar ne 0)] 
;         output[i].ADD_NOISE[0] = add_noise
;     endif

    ;nulls = where(im_conv lt 1.e-13, nc)
    ;if nc gt 0 then ivar[nulls] = 0.
    ;if nc gt 0 then image[nulls] = 0.

    data = {image:image, ivar:ivar, psf:locpsf}
    xcrop=0
    ycrop=0
    imsize = size(data.image, /dimensions)

    makegrid, imsize[0], imsize[1], x, y

    if( not keyword_set(nodvc)) then begin
        diskbulgefit,nod,params,data.image,data.psf,data.ivar,$
          cs,covar,errors,fitstat,dof,bulgeSersic=4.0, $
          nodisk=1, /fixq, _EXTRA={Reff:gals[i].params[1]*1.1, $
                                   q:gals[i].params[3], $
                                   phi:gals[i].params[7]*1.1}
        params[5] += xcrop
        params[6] += ycrop
        output[i].DVCFIT[*,0] = params
        output[i].CHISQ_DVC[0] = cs
        output[i].PERR_DVC[*,0] = errors
        output[i].COVAR_DVC[*,*,0] = covar           
        output[i].MPFIT_STATUS[3,0] = fitstat
        output[i].DOF_DVC[0] = dof
    endif
    
    if( not keyword_set(noexp)) then begin
        diskbulgefit,nod,params,data.image,data.psf,data.ivar,$
          cs,covar,errors,fitstat,dof,bulgeSersic=1.0, $
          nodisk=1, /fixq, _EXTRA={Reff:gals[i].params[1]*1.01, $
                                   q:gals[i].params[3], $
                                   phi:gals[i].params[7]*1.01}
        params[5] += xcrop
        params[6] += ycrop
        output[i].EXPFIT[*,0] = params
        output[i].CHISQ_EXP[0] = cs
        output[i].PERR_EXP[*,0] = errors
        output[i].COVAR_EXP[*,*,0] = covar           
        output[i].MPFIT_STATUS[4,0] = fitstat
        output[i].DOF_EXP[0] = dof
    endif

                                        ;if nosersic is set, then no sersic profile is fit
                                        ; fit a single sersic profile
    if( not keyword_set(nosersic) ) then begin
        diskbulgefit, nod, params, data.image, data.psf, data.ivar, $
          cs, covar, errors, fitstat, dof, $
          /fixq, bulgeSersic=gals[i].params[2]*(1.+0.2*randomn(seed,1)),$;gals[i].params[2]*1.05, $ 
          nodisk=1, freebulge=1, $
          _EXTRA={Reff:gals[i].params[1]*(1.+0.2*randomn(seed,1)), $
                  q:gals[i].params[3], $
                  phi:gals[i].params[7]*1.1}
        
        
        params[5] += xcrop
        params[6] += ycrop
        output[i].SERSICFIT[*,0] = params
        output[i].CHISQ_SERSIC[0] = cs
        output[i].PERR_SERSIC[*,0] = errors
        output[i].COVAR_SERSIC[*,*,0] = covar           
        output[i].MPFIT_STATUS[2,0] = fitstat
        output[i].DOF_SERSIC[0] = dof
    endif


        

endfor
;close,1
;write output fits file: name of file=
name = strsplit(filename, '.', /extract )
file=strsplit(name[0], '/', /extract)
if( n_params() eq 4 ) then begin
    if not keyword_set(add_noise) then add_noise=0.0
    if not keyword_set(add_fwhm) then add_fwhm=0.0
    append = string( start, last-1, add_noise, add_fwhm, $
                     format='("RAWFIT",i02,".",i02,"SN_T",f05.2,"RES",f05.2,".fits")')
endif else append = string( 0, n_elements(gals)-1, format='("RAWFIT",i05,".",i05".fits")')
;outname = name[0] + '/' + file[n_elements(file)-1] + append
outname = outputdir+file[n_elements(file)-1] + append


mwrfits, output, outname, /create
print, 'writing output file: '+outname

END
