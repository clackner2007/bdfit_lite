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
PRO bulge_sn_resolve_sample, filename, start, last, outputdir, $
                       nosersic=nosersic, do_filters=do_filters, redo=redo, $
                       add_fwhm=add_fwhm, add_noise=add_noise
;pref_set, 'idl_cpu_tpool_nthreads', 0, /commit
print, 'cpus: ',!CPU.TPOOL_NTHREADS     

case n_params() of
    2:gals = mrdfits(filename,1)
    4:gals = mrdfits(filename,1,range=[start,last-1])
    else: message, 'Usage: bulgeatlassample, filename of sample'
endcase

;bands array, those with notignore=0 are the ones i'm not going to do
;necessarily the bands in the catalog
bands = ['u','g','r','i','z']
if( keyword_set(do_filters) ) then begin
    notignore = intarr(5)
    for j=0,4 do begin
        hasit = where(do_filters eq bands[j], n)
        if n gt 0 then notignore[j] = 1
    endfor
endif else notignore = [0,0,1,0,0] ; only do r band if nothing is specified
nb = n_elements( where(notignore ne 0) )


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

;insert known values into output strucutre array
output = replicate( output_entry, n_elements(gals) )
struct_assign, gals, output
output.BANDS = reform((where(notignore ne 0))#(intarr(n_elements(gals))+1),n_elements(gals),1)

;file for writing out times
;openw, 1, 'timing.dat'

;loop over each image
for i=0L, n_elements(gals)-1L do begin

    print, gals[i].GALID, ': ', i+1, ' of ', n_elements(gals)
    t1=systime(1)
   ;loop over each band
    for b=0,nb-1 do begin
        
        b_sdss = (where(notignore ne 0))[b]
        print, 'band: ', bands[b_sdss]
 
        
        ;im, read in with 1000 'sfotware bias' and as longs, change to
        ;doubles and subtract bias afterwards
        
        im = sdss_atlas_image(gals[i].RUN, gals[i].CAMCOL,$
                              gals[i].FIELD, gals[i].ID, $
                              rerun=gals[i].RERUN, filter=b_sdss)

        if(n_elements(im) le 1) then begin
            print, 'no atlas image ',gals[i].RUN, gals[i].CAMCOL,$
              gals[i].FIELD, gals[i].ID, gals[i].RERUN
            continue
        endif

        imsize = size(im, /dimensions )
;        print, 'size:',imsize[0],'x',imsize[1]
;         if( imsize[0]*imsize[1] gt 62500 ) then begin
;             print, 'image larger than cutoff size'
;             continue
;         endif
        
        ivar = get_atlas_ivar( im, gals[i].RUN, gals[i].CAMCOL, $
                               gals[i].FIELD, gals[i].ID, $
                               rerun=gals[i].RERUN, filter=b_sdss, $
                               skyCount=sky,gain=gain,sky_err=skyErr,$
                               /use_sky128)

        psffile = string( gals[i].RERUN, gals[i].RUN,$
                          gals[i].CAMCOL, gals[i].RUN, $
                          gals[i].CAMCOL, gals[i].FIELD, $
                          format='("/u/dss/redux/",i0,"/",i0,"/objcs/",i0,"/psField-",i06,"-",i0,"-",i04,".fit")' )
        locpsf = sdss_psf_recon(mrdfits(psffile, b_sdss+1, /silent), $
                                gals[i].ROWC+0.5, $
                                gals[i].COLC+0.5, normalize=1.0 )
          
        im_unbias = (im-1000)*1.0D

        ;convolve the image with a larger psf
        if keyword_set(add_fwhm) then begin
            
            ;blow up the image enough to resolve the psf
            ;resize = ceil(.75/(add_fwhm-1.))*2
            resize=1
            delta_psffwhm = sqrt(add_fwhm^2-1.)*gals[i].PSF_FWHM[b_sdss]

            ;im_blowup = rebin(im_unbias, resize*imsize[0], $
            ;                  resize*imsize[1])
            ;ivar_blowup = rebin(ivar, resize*imsize[0], $
            ;                    resize*imsize[1])
            ;noise_blowup = ivar_blowup
            ;noise_blowup[where(noise_blowup ne 0)] = $            ;  sqrt(1./ivar_blowup[where(noise_blowup ne 0)])
            
           ;psf_blowup = rebin(locpsf, $
           ;                    resize*(size(locpsf,/dimensions))[0], $
           ;                    resize*(size(locpsf, /dimensions))[1])

            ;noise
            noise= ivar
            noise[where(noise ne 0)] = 1./sqrt(noise[where(noise ne 0)])
            
            psf_size=101
            add_psf_correct = PSF_GAUSSIAN(npixel=psf_size, $
                                           fwhm=delta_psffwhm, /normalize)
            im_large = convolve(im_unbias, add_psf_correct $
                                )
            psf_large = convolve(locpsf, add_psf_correct $
                                 )
            noise_large = convolve(noise, add_psf_correct $
                                   )
            
            ;add back noise after smoothing
            ;use the actual noise in the image, ignoring for now the
            ;error in the sky (sky_error) which can be huge
            ; i don't need to add back DC offset errors, but I do need
            ; to add back actual noise
            add_var = (1.-2.*alog(2.)/(!pi*(delta_psffwhm)^2))* $
              (im_unbias+sky)/gain
            posNoise = where(add_var gt 0, pN)
            if pN eq 0 then begin
                print, "could not change noise"
                continue
            endif
            im_large[posNoise] += (randomn(seed, pN))*$
              sqrt(add_var[posNoise])

  
            im_unbias = im_large ;rebin(im_large, imsize[0], imsize[1], /sample)
            im_unbias[where(ivar eq 0)] = 0.
            locpsf = psf_large     ;rebin(psf_large, 51, 51, /sample)
            noise2_large = noise_large^2
            noise2_large[posNoise] += add_var[posNoise]
            ;ivar = ivar
            ivar[where((noise2_large ne 0)and(ivar ne 0))] = $
              1./noise2_large[where((noise2_large ne 0)and(ivar ne 0))]
            ;var = rebin(noise2_large, imsize[0], imsize[1], /sample)
            ;ivarN = dblarr(imsize[0], imsize[1])
            ;ivarN[where(ivar ne 0)] = 1./var[where(ivar ne 0)]      
            output[i].ADD_FWHM[b] = add_fwhm
        endif
        ;add noise to the image, and the ivar
        ;if keyword_set(add_noise) then begin
         ;   seed = systime(1)
          ;  im_unbias[where(ivar ne 0)] += ((randomn(seed, imsize[0], imsize[1]))[where(ivar ne 0)])*$
           ;   sqrt(sky*(add_noise-1.)/gain)
           ; ivar[where(ivar ne 0)] = 1./(1./(ivar[where(ivar ne 0)])+sky*(add_noise-1.)/gain)
           ; output[i].ADD_NOISE[b] = add_noise
        ;endif
        ;increase the total amount noise by upping the sky by add_noise^2
        if keyword_set(add_noise) then begin
            sky_noise = sqrt(sky/gain)
            
            ;seed = systime(1)
            im_unbias[where(ivar ne 0)] += ((randomn(seed, imsize[0], $
                                                     imsize[1]))[where(ivar ne 0)])*$
              sky_noise * $       ;sqrt(1./(ivar[where(ivar ne 0)]))*$
              sqrt(add_noise*add_noise-1.)
            var = ivar
            var[where(ivar ne 0)] = 1./ivar[where(ivar ne 0)] + $
              sky_noise^2*(add_noise^2-1.)
            ivar[where(ivar ne 0)] = 1./var[where(ivar ne 0)] ;((1./(ivar[where(ivar ne 0)]))*(add_noise*add_noise))
            output[i].ADD_NOISE[b] = add_noise
        endif


        data = cropimage(im_unbias, ivar, locpsf, gals[i].PETROR90[b_sdss], $
                         x0=xcrop, y0=ycrop )
        
        ;special case galaxy from reina's sample with bad object subtraction
;         if( gals[i].RUN eq 4263 and gals[i].CAMCOL eq 2 and $
;             gals[i].FIELD eq 69 and gals[i].ID eq 3 ) then begin
;             data2={image:im_unbias[*,10:165], $
;                    ivar:ivar[*,10:165], $
;                    psf:data.psf }
;             data= data2
;             xcrop=0.0
;             ycrop=10.0

;         endif
        imsize = size(data.image, /dimensions)

        ;don't do ridiculously huge images (hack for now)
        ; bad example: SDSSJ091441.16+064610.2
        if( (imsize[0] gt 400 or imsize[1] gt 400) and not $
            (keyword_set(add_fwhm) or keyword_set(add_noise))) then begin
            print, 'too big and image'
            continue
        endif
        
        
;        print, 'cropped size:', imsize[0],'x',imsize[1]
        
         ;data = {image:(im-1000)*1.0D, ivar:ivar, psf:locpsf*1.0D }        
        makegrid, imsize[0], imsize[1], x, y

        ;make the psf image, if it isn't there yet
        if( total(data.psf) lt 0.04 ) then $
          data.psf = makepsf( imsize[0], imsize[1], $
                              gals[i].psf_fwhm[b_sdss] )
        ;normalize psf flux to one
;        data.psfcrop /= total(data.psfcrop)
        if( gals[i].FRACDEV[b_sdss] ge 0.5 ) then $
          reff_guess=gals[i].R_DEV[b_sdss] $
        else reff_guess=gals[i].R_EXP[b_sdss]


        diskbulgefit,nod,params,data.image,data.psf,data.ivar,$
          cs,covar,errors,fitstat,dof,bulgeSersic=4.0, $
          nodisk=1, _EXTRA={Reff:reff_guess, $
                            q:gals[i].AB_DEV[b_sdss], $
                            phi:gals[i].PHI_EXP[b_sdss]}
        params[5] += xcrop
        params[6] += ycrop
        output[i].DVCFIT[*,b] = params
        output[i].CHISQ_DVC[b] = cs
        output[i].PERR_DVC[*,b] = errors
        output[i].COVAR_DVC[*,*,b] = covar           
        output[i].MPFIT_STATUS[3,b] = fitstat
        output[i].DOF_DVC[b] = dof
        
        diskbulgefit,nod,params,data.image,data.psf,data.ivar,$
          cs,covar,errors,fitstat,dof,bulgeSersic=1.0, $
          nodisk=1, _EXTRA={Reff:reff_guess, $
                            q:gals[i].AB_DEV[b_sdss], $
                            phi:gals[i].PHI_EXP[b_sdss]}
        params[5] += xcrop
        params[6] += ycrop
        output[i].EXPFIT[*,b] = params
        output[i].CHISQ_EXP[b] = cs
        output[i].PERR_EXP[*,b] = errors
        output[i].COVAR_EXP[*,*,b] = covar           
        output[i].MPFIT_STATUS[4,b] = fitstat
        output[i].DOF_EXP[b] = dof

        guess = [output[i].EXPFIT[0:4], output[i].EXPFIT[7], $
                 output[i].DVCFIT[0:4], output[i].DVCFIT[7], $
                 output[i].DVCFIT[5:6]-[xcrop,ycrop]]
        guess[6] *= (0.1+gals[i].FRACDEV[b_sdss]*0.8)
        guess[0] *= (1.-(0.1+gals[i].FRACDEV[b_sdss]*0.8))
 
        diskbulgefit, dbparams, bbparams, $
          data.image, data.psf,$
          data.ivar, chsqdb, covardb, errdb, stat, dofdb, bulgeSersic=4.0D, $
          _EXTRA={Reff:reff_guess, $
                  fracdev:gals[i].FRACDEV[b_sdss], $
                  q:gals[i].AB_DEV[b_sdss], $
                  phi:gals[i].PHI_EXP[b_sdss]};,start:guess}
        output[i].MPFIT_STATUS[0,b] = stat

 ;       print, 'timed: ',systime(1)-t1
        guess[8]=1.0
        diskbulgefit, doparams, diparams, $
          data.image, data.psf, $
          data.ivar, chsqdd, covardd, errdd, stat, dofdd,bulgeSersic=1.0D, $
          _EXTRA={Reff:gals[i].R_EXP[b_sdss], $
                  q:gals[i].AB_EXP[b_sdss], $
                  phi:gals[i].PHI_EXP[b_sdss]};,start:guess}
        output[i].MPFIT_STATUS[1,b] = stat


        output[i].DISKFIT[*,b] = [doparams[0:4],doparams[7],$
                                  diparams[0:4],diparams[7],$
                                  diparams[5:6]+[xcrop,ycrop]]
        output[i].BULGEFIT[*,b] = [dbparams[0:4],dbparams[7],$
                                   bbparams[0:4],bbparams[7],$
                                   bbparams[5:6]+[xcrop,ycrop]]
        output[i].PERR_BULGE[*,b] = errdb
        output[i].PERR_DISK[*,b] = errdd
        output[i].COVAR_DISK[*,*,b] = covardd
        output[i].COVAR_BULGE[*,*,b] = covardb
        output[i].CHISQ_BULGE[b] = chsqdb
        output[i].CHISQ_DISK[b] = chsqdd
        output[i].DOF_BULGE[b] = dofdb
        output[i].DOF_DISK[b] = dofdd


        ;if nosersic is set, then no sersic profile is fit
        ; fit a single sersic profile
        if( not keyword_set(nosersic) ) then begin
            diskbulgefit, nod, params, data.image, data.psf, data.ivar, $
              cs, covar, errors, fitstat, dof, bulgeSersic=4.0, $ 
              nodisk=1, freebulge=1, $
              _EXTRA={Reff:reff_guess, $
                      q:gals[i].AB_DEV[b_sdss], $
                      phi:gals[i].PHI_EXP[b_sdss]}
            
            
            params[5] += xcrop
            params[6] += ycrop
            output[i].SERSICFIT[*,b] = params
            output[i].CHISQ_SERSIC[b] = cs
            output[i].PERR_SERSIC[*,b] = errors
            output[i].COVAR_SERSIC[*,*,b] = covar           
            output[i].MPFIT_STATUS[2,b] = fitstat
            output[i].DOF_SERSIC[b] = dof
        endif


        
    endfor

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
