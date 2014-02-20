;+
;____________________________________
;
; Claire Lackner
;  Jan 14 2010
;
; FIT_SAMPLE.PRO
;
; workhorse of fitter, this is what you call to do the fitting
; it writes a file with the output fit parameters. this file will also
;include any column from the input file
;
;ARGUMENTS: 
; filename = input FITS file with list of galaxies you want to fit
;   this list MUST include the following columns:
;   1) NAME unique (string name of object) (for NSA objects use IAU
;name)
;   2) ATLAS_ID (from NSA catalog)
;   3) PARENT_ID (from NSA catalog)
;   These are needed by get_imivarpsf.pro, if you use different
;data/change the data format, you can change these, but you'll
;need to edit get_imivarpsf.pro accordingly
;  4) For the profiles you want to fit (if they aren't default
;ones, see below), you'll need columns PROFILENAME_FIX and
;PROFILENAME_VAL. These contain arrays of 8x(# sersic profiles in fit)
;values and provide the initial conditions (and whether parameters are
;fixed. See make_input.py for details.
;
;start = index in filename of first galaxy to fit (0 -indexed) this
;way you can run multiple instances of the fitter (over condor/pbs)
;for the same input galaxy file and then combine the outputs at the end
;
;last = index in filename table of 1+last galaxy to fit (length of
;file if everything)
;
;outputdir = directory to put output in. If you want to save the model
;images, this needs a subdirectory called models
;
;imagedir = directory where images are. The images have 3 parts, the
;image, the inverse variance plane, and the psf. They go in 3 folders,
;images, ivar, psf
;
;OPTIONS:
;redo = if the input file already includes 'output' columns set this
;option (not often used)
;
;residuals = if you want the code to save images of the fitted models
;galaxies, set this
; 
;crop = if you want the code to crop the images before fitting, set
;this, you won't want this often. It calls cropimage_2.pro.
;
;profiles = list of profiles you want to fit. If you don't set
;this, the code defaults to a sersic profile. If you want the code to
;use the input parameters given in the input file, the names used here
;have to match what's in the file. profiles should be a
;structure, where the keys are the profile names and the values are
;the number of parameters (a multiple of 8), for example,
;profiles={PROFILE1:16, PROFILE2:16}, would fit 2, 2-component
;profiles (using parameters from the input file)
; You can specify default profiles that the code knows and then
;you don't need to set the input parameters, the defaults are:
;DVC, EXP, SER, DVCEXP
;
;filter='u','g','r',...name of the band you are using (SDSS bands for
;NSA)
;
;freesky = set this if you want the sky to be a free parameter
;
;cutoff = set this if you want the models to smoothly go to zero after
;4 scalelengths (exponential) 8 scalelengths (dvc or sersic). This is
;the typical SDSS behavior. It will slightly change results, so be
;consistent
;
;debug = set this if you want the code to print extra information
;
;OUTPUTS: the program writes a file called
;outputdir/RAWFITxxx{start}.xxx{last}.fits this file has all the
;columns from the input file plus the following:
;XCROP, YCROP  - position of left, bottom corner if cropping is used,
;                otherwise, 0,0
;XLEN, YLEN - size of image fit (smaller than original if cropping)
;For each profile fit there are the following columns
;FIT_PROFILENAME: final fit parameters for each Sersic component fit
;there are 8 parameters (this goes for the input to) they are in this
;order (see sersic.pro for formulas)
;1) surface brightness at half-light radius
;2) half-light radius (refers to semi-major axis)
;3) sersic index
;4) axis ratio of profile (minor/major)
;5) boxiness/diskiness by default I set this to zero, so should you
;6) the x coordinate of the central position of the profiles, by
;default all components have the same center
;7) the y coordinate of the central position
;8) the position angle of the component (different components can be
;at different angles) in radians (counterclockwise from the x axis)
;If the profile to be fit has 2 components, FIT_XXX is a 16-element
;array for each galaxy in the list, 
;CHISQ_PROFILENAME: reduced chi^2 of fit
;COVAR_PROFILENAME: covariance matrix of fit parameters (size of FIT^2)
;PERR_PROFILENAME: errors (1-sigma) in fit parameters (same size as FIT_XX)
;STAT_PROFILENAME: status of fit, anything other than 0 or 5 is ok
;DOF_PROFILENAME: degrees of freedom in fit (#pixels - #parameters in fit)
;SKY_PROFILENAME: sky value (0 if sky isn't free)
;SKYERR_PROFILENAME: error in sky value (note there are also
;covariances one could return
;FLUXRATIO_PROFILENAME: flux ratio of first (usually highest sersic
;index/smallest component) in a profile to the total flux in the
;profile (bulge to total ratio)
;
;if you turn residuals on, then the program also saves a copy of the
;best-fit model. It does this to a file outputdir/models/M{NAME}.fits
;One file for each galaxy fit. So the profiles in there are in the
;order you fit them (you'll have to pay attention)
;____________________________________
;-

PRO fit_sample, filename, start, last, outputdir, imagedir, $
                redo=redo, residuals=residuals, $
                crop=crop, profiles=profiles, $
                filter=filter, $
                freesky=freesky, cutoff=cutoff, debug=debug
print, 'cpus: ',!CPU.TPOOL_NTHREADS     

gals = mrdfits(filename,1,range=[start,last-1])

if (n_elements(profiles) eq 0) then begin
    message, 'ERROR: No profiles to fit'
    ;default is a single sersic profile
    ;profiles={SER:8}
endif
;profiles = hash(profiles, /extract)

;default filter is 'r'
if (n_elements(filter) eq 0) then filter='r'
filters=['u', 'g', 'r', 'i', 'z']
if (where(filters eq filter) eq -1) then filter='r'
band=(where(filters eq filter))[0]

if not keyword_set(redo) then begin

   output_entry = create_struct(gals[0], ['XCROP', 'YCROP',$
                                           'XLEN', 'YLEN'], $
                               0L, 0L, 0L, 0L)

   for indp=0, n_tags(profiles)-1 do begin
      prof = (tag_names(profiles))[indp]
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
                   'SKY_'+prof, 'SKYERR_'+prof, 'FLUXRATIO_'+prof]
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

    print, gals[i].NAME, ': ', i+1, ' of ', n_elements(gals)
    t1=systime(1)
    
    ;get the image, ivar, and psf
    data = get_imivarpsf(gals[i].NAME, imagedir, $
                         gals[i].ATLAS_ID, gals[i].PARENT_ID, $
                         band=band)

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

    ;if residuals is set, then save the model images
    if keyword_set(residuals) then begin
       modName = string('models/M'+output[i].NAME+'.fits')
       modName = outputdir + modName
    endif

    ;do the profile fitting
    for ip=0, n_tags(profiles)-1 do begin
       prof = (tag_names(profiles))[ip]
       nvar = profiles.(ip)
       plist = indgen(nvar/8)*8
       if keyword_set(debug) then print, "doing profile "+prof

                                ;check if we allow a freesky otherwise
                                ;it defaults to 0.0 TODO: allow the
                                ;sky to be fixed to a value other than
                                ;0, and let this change for different
                                ;profiles 
       if keyword_set(freesky) then sky=0.0D else sky=!NULL

                                ;check to see if there are fixed
                                ;parameters set in the input the names
                                ;of fixed parameters should be
                                ;{PROF}_FIX which is a boolean array
                                ;(=1 for fixed) and {PROF}_VAL which
                                ;is any array of starting values if
                                ;neither are set, check and see what
                                ;the profile name is, if its 'DVC',
                                ;'EXP', 'SER', 'DVCEXP' then do the
                                ;standard things
       tname1 = where(strcmp(outnames, prof+'_FIX') eq 1)
       tname2 = where(strcmp(outnames, prof+'_VAL') eq 1)

       if ((tname1 ne -1) and (tname2 ne -1)) then begin
          fixed_params = output[i].(tname1)
          tempparams = output[i].(tname2)
          tempparams[plist+5] -= xcrop
          tempparams[plist+6] -= ycrop
          ;set the initial scalings and central position
          params = init_cond(tempparams, fixed_params, data.image)
          if keyword_set(debug) then print, "using giving IC"
       endif else begin
          case prof of
             'DVC': begin
                 temp = default_init_cond(1,data.image, sersics=[4.0])
                 params = temp.params
                 fixed_params = temp.fixed
                 fixed_params[2+plist] = 1
                 if keyword_set(debug) then print, "using default DVC IC"
              end
             'EXP': begin
                temp = default_init_cond(1,data.image, sersics=[1.0])
                params = temp.params
                fixed_params = temp.fixed
                fixed_params[2+plist] = 1
                if keyword_set(debug) then print, "using default EXP IC"
             end
             'SER': begin
                 temp = default_init_cond(1,data.image, sersics=[4.0])
                 params = temp.params
                 fixed_params = temp.fixed
                 if keyword_set(debug) then print, "using default SER IC"
              end
             'DVCEXP': begin
                temp = default_init_cond(2, data.image, sersic=[4.0,1.0])
                params = temp.params
                fixed_params = temp.fixed
                fixed_params[2+plist] = 1
                if keyword_set(debug) then print, "using default DVCEXP IC"
             end
             else: begin
                temp = default_init_cond(nvar/8, data.image)
                params = temp.params
                fixed_params = temp.fixed
                if keyword_set(debug) then print, "using default IC"
             end
          endcase
       endelse 

       db_flexfit, params, data.image, data.psf, data.ivar, $
                   chsq, covar, err, stat, dof, skyVal=sky, $
                   free_sky=freesky, fix_params=fixed_params, $
                   debug=debug

       params[plist+5] += xcrop
       params[plist+6] += ycrop

       output[i].(WHERE(STRCMP(outnames,'FIT_'+prof) EQ 1)) = params
       output[i].(WHERE(STRCMP(outnames,'COVAR_'+prof) EQ 1)) = $
          covar[0:nvar-1, 0:nvar-1]
       output[i].(WHERE(STRCMP(outnames,'PERR_'+prof) EQ 1)) = err[0:nvar-1]
       output[i].(WHERE(STRCMP(outnames,'CHISQ_'+prof) EQ 1)) = chsq
       output[i].(WHERE(STRCMP(outnames,'STAT_'+prof) EQ 1)) = stat
       output[i].(WHERE(STRCMP(outnames,'DOF_'+prof) EQ 1)) = dof
       output[i].(WHERE(STRCMP(outnames,'SKY_'+prof) EQ 1)) = sky
       output[i].(WHERE(STRCMP(outnames,'SKYERR_'+prof) EQ 1)) = err[nvar]
       if n_elements(params) gt 8 then $
          output[i].(WHERE(STRCMP(outnames,'FLUXRATIO_'+prof) EQ 1)) = $
          bulgetotot(params, cutoff=cutoff) $
       else $
          output[i].(WHERE(STRCMP(outnames,'FLUXRATIO_'+prof) EQ 1)) = 1.0

       if keyword_set(residuals) then begin
          model = model_image(origsize[0], origsize[1], params, $
                              data.psf, cutoff=cutoff)
          mwrfits, model, modName
       endif
       
    endfor
   
 endfor
;close,1
;write output fits file: name of file=
append = string( start, last-1, format='("RAWFIT",i05,".",i05,".fits")')
outname = outputdir+append

mwrfits, output, outname, /create
print, 'writing output file: '+outname

END
