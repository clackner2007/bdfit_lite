;+
;____________________________________
;
; Claire Lackner
;  Jan 14 2010
;
; FIT_SAMPLE.PRO
;
; fits a set of profiles to a set of galaxies  
;
; 
;____________________________________
;-

PRO fit_sample, filename, start, last, outputdir, imagedir, $
                redo=redo, residuals=residuals, $
                crop=crop, profiles=profiles, $
                freesky=freesky, cutoff=cutoff
print, 'cpus: ',!CPU.TPOOL_NTHREADS     

gals = mrdfits(filename,1,range=[start,last-1])

;default is a single sersic profile
if (n_elements(profiles) eq 0) then profiles={SER:8}
;profiles = hash(profiles, /extract)


if not keyword_set(redo) then begin

   output_entry = create_struct(gals[0], ['XCROP', 'YCROP',$
                                           'XLEN', 'YLEN'], $
                               0L, 0L, 0L, 0L)

   for indp=0, n_tags(profiles) do begin
      prof = (get_tags(profiles))[indp]
      nfit = profiles.(indp)
                                ;output entries for each profile
                                ;include: the fit parameters, the
                                ;reduced chi^2 value, the 1-sigma
                                ;fitting errors in the parameters, the
                                ;covariance matrix of the parameters,
                                ;the status of the fitter, the degrees
                                ;of freedom, the sky value, the error
                                ;in the sky, and the flux ratio of the
                                ;profiles fit. This is ALWAYS
                                ;first/total (like B/T)
      ent_names = ['FIT_'+prof, 'CHISQ_'+prof, $
                   'COVAR_'+prof, 'PERR_'+prof, $
                   'STAT_'+prof, 'DOF_'+prof, $
                   'SKY_'+prof, 'SKYERR_'+prof, 'FLUXRAITO_'+prof]
      output_entry = create_struct(output_entry, $
                                   ent_names, $
                                   dblarr(nfit), 0.0D, $
                                   dblarr(nfit,nfit), dblarr(nfit), $
                                   0, 0L, 0.0D, 0.0D, 0.0D)                             
                                   
   endfor

endif else output_entry = create_struct(gals[0])

;insert known values into output strucutre array
output = replicate( output_entry, n_elements(gals) )
struct_assign, gals, output
outnames = tag_names(output)


;loop over each image
for i=0L, n_elements(gals)-1L do begin

    print, gals[i].IDENT, ': ', i+1, ' of ', n_elements(gals)
    t1=systime(1)
    
    ;get the image, ivar, and psf
    data = get_imivarpsf(gals[i].IDENT, strtrim(gals[i].FILENAME), imagedir)

                                ;figure out the size of the image and
                                ;do cropping (if wanted, you rarely
                                ;want this)
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

    ;save the image size to the structure (useful, esp. if cropped)
    output[i].XLEN=imsize[0]
    output[i].YLEN=imsize[1]
    output[i].XCROP=xcrop
    output[i].YCROP=ycrop

            

    ;do the profile fitting
    for ip=0, n_tags(profiles) do begin
       prof = (tag_names(profiles))[ip]
       nvar = profiles.(ip)
       plist = indgen(nvar/8)*8

                                ;check if we allow a freesky otherwise
                                ;it defaults to 0.0 TODO: allow the
                                ;sky to be fixed to a value other than
                                ;0, and let this change for different
                                ;profiles 
       if keyword_set(freesky) then sky=0.0D else sky=!NULL 


       db_flexfit, params, data.image, data.psf, data.ivar, $
                   chsq, covar, err, stat, dof, skyVal=sky, $
                   free_sky=freesky
       params[plist+5] += xcrop
       params[plist+6] += ycrop

       output[i].(WHERE(STRCMP(outnames,'FIT_'+prof) EQ 1)) = params
       output[i].(WHERE(STRCMP(outnames,'COVAR_'+prof) EQ 1)) = covar[0:nvar-1, 0:nvar-1]
       output[i].(WHERE(STRCMP(outnames,'PERR_'+prof) EQ 1)) = err[0:nvar-1]
       output[i].(WHERE(STRCMP(outnames,'CHISQ_'+prof) EQ 1)) = chsq
       output[i].(WHERE(STRCMP(outnames,'STAT_'+prof) EQ 1)) = stat
       output[i].(WHERE(STRCMP(outnames,'DOF_'+prof) EQ 1)) = dof
       output[i].(WHERE(STRCMP(outnames,'SKY_'+prof) EQ 1)) = sky
       output[i].(WHERE(STRCMP(outnames,'SKYERR_'+prof) EQ 1)) = err[nvar]
       output[i].(WHERE(STRCMP(outnames,'FLUXRATIO_'+prof) EQ 1)) = $
          bulgetotot(params, cutoff=cutoff)

    endfor
   
endfor
;close,1
;write output fits file: name of file=
append = string( start, last-1, format='("RAWFIT",i05,".",i05,".fits")')
outname = outputdir+append

mwrfits, output, outname, /create
print, 'writing output file: '+outname

END
