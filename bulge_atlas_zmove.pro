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
PRO bulge_atlas_zmove, filename, start, last, outputdir, imfolder, $
                       nosersic=nosersic, do_filters=do_filters, $
                       change_res=change_res, redo=redo, gausspsf=gausspsf
print, 'cpus: ',!CPU.TPOOL_NTHREADS     

case n_params() of
    3:gals = mrdfits(filename,1)
    5:gals = mrdfits(filename,1,range=[start,last-1])
    else: message, 'Usage: bulgeatlas_zmove, filename of sample'
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
    output_entry = create_struct( gals[0], ['CURRENT_Z', 'CURR_FLUX_RATIO',$
                                            'CURR_SIZE_RATIO',$
                                            'BULGEFIT', 'DISKFIT',$
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
                                            'BANDS', 'LORES_BULGEFIT', $
                                            'LORES_DISKFIT', 'LORES_SERSIC',$
                                            'LORES_DVC', 'LORES_EXP'], $
                                  0.0, 0.0, 0.0, $
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
                                  intarr(nb), dblarr(14,nb), $
                                  dblarr(14,nb), dblarr(8,nb),$
                                  dblarr(8,nb), dblarr(8,nb) )
endif else output_entry = create_struct(gals[0],  ['CURRENT_Z', 'CURR_FLUX_RATIO',$
                                                   'CURR_SIZE_RATIO',$
                                                   'LORES_BULGEFIT', $
                                                   'LORES_DISKFIT', 'LORES_SERSIC',$
                                                   'LORES_DVC', 'LORES_EXP'], $
                                        0.0, 0.0, 0.0, $
                                        dblarr(14,nb), $
                                        dblarr(14,nb), dblarr(8,nb),$
                                        dblarr(8,nb), dblarr(8,nb) )

;insert known values into output strucutre array
nreds = n_elements(gals[0].MOVED_Z)
output = replicate( output_entry, n_elements(gals)*(nreds) )

;file for writing out times
;openw, 1, 'timing.dat'

;loop over each image
for i=0L, n_elements(gals)-1L do begin

    print, gals[i].GALID, ': ', i+1, ' of ', n_elements(gals)
    t1=systime(1)

;loop over each redshift step
    for red=0L, nreds-1 do begin
        if gals[i].MOVED_Z[red] lt 1.0e-6 then continue

        print, 'Z: ',gals[i].MOVED_Z[red]
        struct_assign, gals[i], output_entry
        output[nreds*i+red] = output_entry
        output[nreds*i+red].CURRENT_Z = gals[i].MOVED_Z[red]
        output[nreds*i+red].CURR_FLUX_RATIO = gals[i].FLUX_RATIO[red]
            
        size_ratio = lumdist(gals[i].Z,/silent)/lumdist(gals[i].MOVED_Z[red],/silent)*$
          (1.0+gals[i].MOVED_Z[red])^2/(1.0+gals[i].Z)^2
        
        output[nreds*i+red].CURR_SIZE_RATIO = size_ratio
            
   ;loop over each band
    for b=0,nb-1 do begin
        
        b_sdss = (where(notignore ne 0))[b]
        print, 'band: ', bands[b_sdss]
 
        
        ;im, read in with 1000 'sfotware bias' and as longs, change to
        ;doubles and subtract bias afterwards
        path=imfolder+gals[i].GALID+"-Z-"

        im_file=path+string(red,format='(i02,".fits")')

        im = mrdfits(im_file,1,/silent)
        ivar = mrdfits(im_file,2,/silent)

        imsize = size(im, /dimensions )
        if n_elements(im) lt 4 then continue
;        print, 'size:',imsize[0],'x',imsize[1]
;         if( imsize[0]*imsize[1] gt 62500 ) then begin
;             print, 'image larger than cutoff size'
;             continue
;         endif


        ;if not noted to use the gaussain psf, use the standard sdss psf
        if not keyword_set(gausspsf) then begin
            
            psffile = string( gals[i].RERUN, gals[i].RUN,$
                              gals[i].CAMCOL, gals[i].RUN, $
                              gals[i].CAMCOL, gals[i].FIELD, $
                              format='("/u/dss/redux/",i0,"/",i0,"/objcs/",i0,"/psField-",i06,"-",i0,"-",i04,".fit")' )
            locpsf = sdss_psf_recon(mrdfits(psffile, b_sdss+1, /silent), $
                                    gals[i].ROWC+0.5, $
                                    gals[i].COLC+0.5, normalize=1.0 )
        endif else begin
            ;use a gaussian for the psf for fake images
            locpsf = PSF_GAUSSIAN(npixel=51,$
                                  fwhm=gals[i].PSF_FWHM[b_sdss], $
                                  /normalize)
        endelse
        
            
        im_unbias = (im-1000)*1.0D

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
        if( imsize[0] gt 400 or imsize[1] gt 400 ) then begin
            print, 'too big and image'
            continue
        endif
         
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
 
        diskbulgefit, dbparams, bbparams, $
          data.image, data.psf,$
          data.ivar, chsqdb, covardb, errdb, stat, dofdb, bulgeSersic=4.0D, $
          _EXTRA={Reff:reff_guess*size_ratio, $
                  fracdev:gals[i].FRACDEV[b_sdss], $
                  q:gals[i].AB_DEV[b_sdss], $
                  phi:gals[i].PHI_EXP[b_sdss]}
        output[nreds*i+red].MPFIT_STATUS[0,b] = stat

 ;       print, 'timed: ',systime(1)-t1
        diskbulgefit, doparams, diparams, $
          data.image, data.psf, $
          data.ivar, chsqdd, covardd, errdd, stat, dofdd,bulgeSersic=1.0D, $
          _EXTRA={Reff:gals[i].R_EXP[b_sdss]*size_ratio, $
                  fracdev:0.0, $
                  q:gals[i].AB_EXP[b_sdss], $
                  phi:gals[i].PHI_EXP[b_sdss]}
        output[nreds*i+red].MPFIT_STATUS[1,b] = stat


        output[nreds*i+red].DISKFIT[*,b] = [doparams[0:4],doparams[7],$
                                  diparams[0:4],diparams[7],$
                                  diparams[5:6]+[xcrop,ycrop]]
        output[nreds*i+red].BULGEFIT[*,b] = [dbparams[0:4],dbparams[7],$
                                   bbparams[0:4],bbparams[7],$
                                   bbparams[5:6]+[xcrop,ycrop]]
        output[nreds*i+red].PERR_BULGE[*,b] = errdb
        output[nreds*i+red].PERR_DISK[*,b] = errdd
        output[nreds*i+red].COVAR_DISK[*,*,b] = covardd
        output[nreds*i+red].COVAR_BULGE[*,*,b] = covardb
        output[nreds*i+red].CHISQ_BULGE[b] = chsqdb
        output[nreds*i+red].CHISQ_DISK[b] = chsqdd
        output[nreds*i+red].DOF_BULGE[b] = dofdb
        output[nreds*i+red].DOF_DISK[b] = dofdd


        ;if nosersic is set, then no sersic profile is fit
        ; fit a single sersic profile
        if( not keyword_set(nosersic) ) then begin
            diskbulgefit, nod, params, data.image, data.psf, data.ivar, $
              cs, covar, errors, fitstat, dof, bulgeSersic=4.0, $ 
              nodisk=1, freebulge=1, $
              _EXTRA={Reff:reff_guess*size_ratio, $
                      q:gals[i].AB_DEV[b_sdss], $
                      phi:gals[i].PHI_EXP[b_sdss]}
            
            
            params[5] += xcrop
            params[6] += ycrop
            output[nreds*i+red].SERSICFIT[*,b] = params
            output[nreds*i+red].CHISQ_SERSIC[b] = cs
            output[nreds*i+red].PERR_SERSIC[*,b] = errors
            output[nreds*i+red].COVAR_SERSIC[*,*,b] = covar           
            output[nreds*i+red].MPFIT_STATUS[2,b] = fitstat
            output[nreds*i+red].DOF_SERSIC[b] = dof
        endif


        diskbulgefit,nod,params,data.image,data.psf,data.ivar,$
          cs,covar,errors,fitstat,dof,bulgeSersic=4.0, $
          nodisk=1, _EXTRA={Reff:reff_guess*size_ratio, $
                            q:gals[i].AB_DEV[b_sdss], $
                            phi:gals[i].PHI_EXP[b_sdss]}
        params[5] += xcrop
        params[6] += ycrop
        output[nreds*i+red].DVCFIT[*,b] = params
        output[nreds*i+red].CHISQ_DVC[b] = cs
        output[nreds*i+red].PERR_DVC[*,b] = errors
        output[nreds*i+red].COVAR_DVC[*,*,b] = covar           
        output[nreds*i+red].MPFIT_STATUS[3,b] = fitstat
        output[nreds*i+red].DOF_DVC[b] = dof
        
        diskbulgefit,nod,params,data.image,data.psf,data.ivar,$
          cs,covar,errors,fitstat,dof,bulgeSersic=1.0, $
          nodisk=1, _EXTRA={Reff:reff_guess*size_ratio, $
                            q:gals[i].AB_DEV[b_sdss], $
                            phi:gals[i].PHI_EXP[b_sdss]}
        params[5] += xcrop
        params[6] += ycrop
        output[nreds*i+red].EXPFIT[*,b] = params
        output[nreds*i+red].CHISQ_EXP[b] = cs
        output[nreds*i+red].PERR_EXP[*,b] = errors
        output[nreds*i+red].COVAR_EXP[*,*,b] = covar           
        output[nreds*i+red].MPFIT_STATUS[4,b] = fitstat
        output[nreds*i+red].DOF_EXP[b] = dof


        if( keyword_set(change_res) ) then begin
            orig_data = get_imivarpsf(gals[i].RUN, gals[i].CAMCOL,$
                                      gals[i].FIELD, gals[i].ID, $
                                      gals[i].COLC[2], $
                                      gals[i].ROWC[2], $
                                      rerun=gals[i].RERUN )
            ;fhwm of gaussian to convolve with to lower the resolution
            ;to the resolution of the higher z image

            osize = size(orig_data.image, /dimensions)
            low_res_psf = makepsf(osize[0], osize[1], 1./size_ratio)
            makegrid, osize[0], osize[1], x, y

            orig_data.image = convolve(orig_data.image, low_res_psf,$
                                       ft_psf=ftpsf)
            orig_data.ivar = convolve(orig_data.ivar, low_res_psf,$
                                      ft_psf=ftpsf)
            
            diskbulgefit, dbparams, bbparams, $
              orig_data.image, orig_data.psf,$
              orig_data.ivar, chsqdb, covardb, errdb, stat, dofdb, bulgeSersic=4.0D, $
              _EXTRA={Reff:reff_guess, $
                      fracdev:gals[i].FRACDEV[b_sdss], $
                      q:gals[i].AB_DEV[b_sdss], $
                      phi:gals[i].PHI_EXP[b_sdss]}
            diskbulgefit, doparams, diparams, $
              orig_data.image, orig_data.psf,$
              orig_data.ivar, chsqdb, covardb, errdb, stat, dofdb, bulgeSersic=1.0D, $
              _EXTRA={Reff:reff_guess, $
                      fracdev:gals[i].FRACDEV[b_sdss], $
                      q:gals[i].AB_DEV[b_sdss], $
                      phi:gals[i].PHI_EXP[b_sdss]}
            output[nreds*i+red].LORES_BULGEFIT[*,b] = [dbparams[0:4],dbparams[7],$
                                             bbparams[0:4],bbparams[7],$
                                             bbparams[5:6]]
            output[nreds*i+red].LORES_DISKFIT[*,b] = [doparams[0:4],doparams[7],$
                                                      diparams[0:4],diparams[7],$
                                                      diparams[5:6]+[xcrop,ycrop]]
            diskbulgefit, nod, params, orig_data.image, orig_data.psf, orig_data.ivar, $
              cs, covar, errors, fitstat, dof, bulgeSersic=4.0, $ 
              nodisk=1, freebulge=1, $
              _EXTRA={Reff:reff_guess, $
                      q:gals[i].AB_DEV[b_sdss], $
                      phi:gals[i].PHI_EXP[b_sdss]}
            
             output[nreds*i+red].LORES_SERSIC[*,b] = params

             diskbulgefit,nod,params,orig_data.image,orig_data.psf,orig_data.ivar,$
               cs,covar,errors,fitstat,dof,bulgeSersic=4.0, $
               nodisk=1, _EXTRA={Reff:reff_guess, $
                                 q:gals[i].AB_DEV[b_sdss], $
                                 phi:gals[i].PHI_EXP[b_sdss]}
             output[nreds*i+red].LORES_DVC[*,b] = params
        
             diskbulgefit,nod,params,orig_data.image,orig_data.psf,orig_data.ivar,$
               cs,covar,errors,fitstat,dof,bulgeSersic=1.0, $
               nodisk=1, _EXTRA={Reff:reff_guess, $
                                 q:gals[i].AB_DEV[b_sdss], $
                                 phi:gals[i].PHI_EXP[b_sdss]}
             output[nreds*i+red].LORES_EXP[*,b] = params
             
        endif
       
        
    endfor
endfor 

endfor

;output.BANDS = (where(notignore ne 0))#(intarr(n_elements(output))+1)
output.BANDS = (intarr(n_elements(output))+1)#(where(notignore ne 0))

;close,1
;write output fits file: name of file=
name = strsplit(filename, '.', /extract )
file=strsplit(name[0], '/', /extract)
if( n_params() eq 5 ) then $
  append = string( start, last-1, format='("RAWFIT",i05,".",i05,".fits")') $
else append = string( 0, n_elements(gals)-1, format='("RAWFIT",i05,".",i05,".fits")')
;outname = name[0] + '/' + file[n_elements(file)-1] + append
outname = outputdir+file[n_elements(file)-1] + append


mwrfits, output, outname, /create
print, 'writing output file: '+outname

END
