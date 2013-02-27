;____________________________________
;
; Claire Lackner
;  Jan 14 2010
;
;  bulgefit.pro
;
;  attempts to fit a bulge+disk profile to the VAGC lowz catalog data
;in a way similar to testbdecomp does for artificially generated images
;
; 
;____________________________________


;uses gals which is list of galaxies from sdss lowz catalog
;takes 4 arguments, the ranges in which to look for sources in the
;lowz catalog, both given in degrees
PRO bulgefit, ramin, ramax, decmin, decmax, fit_size

;usage statement
if( n_params() ne 5 ) then message, $
  'Usage: sersicfit RA_min RA_max DEC_min DEC_max (degrees) [array of circles to fit]'

;data is (image, noise, psf) made public so function can access psf
;common datastruct, data

;read in the galaxies
gals = lowz([ramin,ramax], [decmin,decmax], status=s)

if( s eq -1 ) then message, 'no galaxies found'

;bands array, those with notignore=0 are the ones i'm not going to do
;necessarily the bands in the catalog
bands = ['u','g','r','i','z']
notignore = [0,1,1,1,0]
nb = n_elements( bands )

;number of different fitsizes, excluding entire galaxy fit
;each fitsize, given in the header is a fraction of the Petrosian
;radius
nfits=n_elements(fit_size)



;output entry for each object, the size of the array is for the bands,
;the default is five
output_entry = { NAME:gals.name[0], RA:gals.ra[0], DEC:gals.dec[0], $
                 RUN:gals.run[0], CAMCOL:gals.camcol[0], $
                 FIELD:gals.field[0], ID:gals.id[0], $
                 RERUN:gals.rerun[0], $
;                  DOF:dblarr(nb), $
                 CHI_SQUARE_BULGE:dblarr(nb), $
                 CHI_SQUARE_DISK:dblarr(nb), $
                 DISKFIT:dblarr(8,nb), $
                 BULGEFIT:dblarr(8,nb), $
                 FITRANGE:fit_size, NORM:dblarr(nfits+1,nb), $
                 SERSIC_N:dblarr(nfits+1,nb), $
                 R50:dblarr(nfits+1,nb), $
                 ELLIP:dblarr(nfits+1,nb), C:dblarr(nfits+1,nb), $
                 XCEN:dblarr(nfits+1,nb), YCEN:dblarr(nfits+1,nb), $
                 PHI:dblarr(nfits+1,nb), $
                 CHI_SQUARE:dblarr(nfits+1,nb), $
;                  CHI_SQUARE:dblarr(nfits+1), DOF:dblarr(nfits+1), $
;                  COVAR:dblarr(64,5), $
;                  SERSIC_NORM:dblarr(nb), SERSIC_NORM_ERR:dblarr(nb), $
;                  SERSIC_N:dblarr(nb), SERSIC_N_ERR:dblarr(nb), $
;                  SERSIC_R50:dblarr(nb), SERSIC_R50_ERR:dblarr(nb), $
;                  SERSIC_FLUX:dblarr(nb), SERSIC_FLUX_ERR:dblarr(nb), $
;                  B_A:dblarr(nb), B_A_ERR:dblarr(nb), $
;                  C_ELLIP:dblarr(nb), C_ELLIP_ERR:dblarr(nb), $
;                  X_CENTER:dblarr(nb), X_CENTER_ERR:dblarr(nb), $
;                  Y_CENTER:dblarr(nb), Y_CENTER_ERR:dblarr(nb), $
;                  PHI_MAJ:dblarr(nb), PHI_MAJ_ERR:dblarr(nb), $
                 BLANTON_N:dblarr(nb), BLANTON_FLUX:dblarr(nb), $
                 BLANTON_R50:dblarr(nb), SDSS_B_A:0.0, $
                 PETRORAD:0.0D $
;                  FITSIZE:dblarr(nfits), FS_SERSIC_N:dblarr(nfits, nb), $
;                  FS_SERSIC_R50:dblarr(nfits, nb), FS_B_A:dblarr(nfits, nb) $
               }

;insert known values into output strucutre array
output = replicate( output_entry, n_elements(gals.name) )

output[*].NAME = gals.name
output[*].RA = gals.ra
output[*].DEC = gals.dec
output[*].RUN = gals.run
output[*].CAMCOL = gals.camcol
output[*].FIELD = gals.field
output[*].ID = gals.id
output[*].RERUN = gals.rerun
output[*].BLANTON_N = gals.sersic_n
output[*].BLANTON_FLUX = gals.flux
output[*].BLANTON_R50 = gals.reff*0.396
output[*].SDSS_B_A = gals.ba_sdss
output[*].PETRORAD = gals.RPETRO

;loop over each image
for i=0L, n_elements(gals.name)-1L do begin

    print, gals.name[i]
    print, i+1, ' of ', n_elements(gals.name)
    
   ;loop over each band
    for b=0,n_elements(bands)-1 do begin
        
        if( notignore[b] eq 0) then begin
;             output[i].DOF[b] = 0
;             output[i].CHI_SQUARE[b] = 0
;             output[i].SERSIC_NORM[b] = 0
;             output[i].SERSIC_NORM_ERR[b] = 0
;             output[i].SERSIC_N[b] = 0
;             output[i].SERSIC_N_ERR[b] = 0
;             output[i].SERSIC_R50[b] = 0
;             output[i].SERSIC_R50_ERR[b] = 0
;             output[i].B_A[b] = 0
;             output[i].B_A_ERR[b] = 0
;             output[i].C_ELLIP[b] = 0
;             output[i].C_ELLIP_ERR[b] = 0
;             output[i].X_CENTER[b] = 0
;                                 ;in pixels with center indexed by integer
;             output[i].X_CENTER_ERR[b] = 0
;             output[i].Y_CENTER[b] = 0
;             output[i].Y_CENTER_ERR[b] = 0
;                                 ;phi in radians
;             output[i].PHI_MAJ[b] = 0
;             output[i].PHI_MAJ_ERR[b] = 0
;             output[i].COVAR[*,b] = 0
            
;             output[i].SERSIC_FLUX[b] = 0
            continue
        endif

        ;structure containing, image and inverse variance and psf, 
        data = readimages( gals.folder[i], gals.name[i], bands[b], $
                         status=status)
        
        ;check to make sure the file was properly read
        if( status lt 0 ) then continue
        
        ;crop unnecessary noise on outside of image
        newdata = cropimage( data.image, data.ivar, (gals.rpetro)[i], $
                             x0=xcrop, y0=ycrop )
  
        ;combine new and old data
        data = create_struct( data, newdata )

        imsize = size(data.imagecrop, /dimensions)
        makegrid, imsize[0], imsize[1], x, y

        ;make the psf image, if it isn't there yet
        if( total(data.psfcrop) lt 0.04 ) then $
          data.psfcrop = makepsf( imsize[0], imsize[1], $
                                  (gals.psf_fwhm)[b,i] )
        ;normalize psf flux to one
        data.psfcrop /= total(data.psfcrop)
 
        diskbulgefit, dbparams, bbparams, $
          data.imagecrop, (gals.psf_fwhm)[b,i],$
          data.ivarcrop, chsqdb, bulgeSersic=4.0, $
          _EXTRA={Reff:(gals.reff)[b,i]*0.396, $
                  q:(gals.ba_sdss)[i]}
        diskbulgefit, doparams, diparams, $
          data.imagecrop, (gals.psf_fwhm)[b,i],$
          data.ivarcrop, chsqdd, bulgeSersic=1.0, $
          _EXTRA={Reff:(gals.reff)[b,i]*0.396, $
                  q:(gals.ba_sdss)[i]}

        if( chsqdb lt chsqdd ) then begin
            dparams = dbparams
            bparams = bbparams
        endif else begin
            dparams = doparams
            bparams = diparams
        endelse
        
        if( dparams[0]/bparams[0] lt 1.0e-2 or $
            dparams[1] lt 1.0e-1 or $
            bparams[3] lt 0.05 ) then $
          subtractImage = data.imagecrop $
        else begin $
          subtractImage = data.imagecrop - pixelfluxpsf( x, y,dparams, $
                                                         _EXTRA={psfImage:data.psfcrop})
            print, 'using disk-subtracted image'
        endelse
;        subtractImage = data.imagecrop
;         fitImage, val, fit_size, subtractImage, $
;           (gals.psf_fwhm)[b,i], data.ivarcrop
        

        output[i].DISKFIT[*,b] = dparams
        output[i].BULGEFIT[*,b] = bparams
        output[i].CHI_SQUARE_BULGE[b] = chsqdb
        output[i].CHI_SQUARE_DISK[b] = chsqdd
;         output[i].NORM[*,b] = val[*,0]
;         output[i].SERSIC_N[*,b] = val[*,2]
;         output[i].R50[*,b] = val[*,1]
;         output[i].ELLIP[*,b] = val[*,3]
;         output[i].C[*,b] = val[*,4]
;         output[i].XCEN[*,b] = val[*,5]
;         output[i].YCEN[*,b] = val[*,6]
;         output[i].PHI[*,b] = val[*,7]
;         output[i].CHI_SQUARE[*,b] = val[*,9]/val[*,10]
        
;         center model at highest value pixel
;         indmax = where( data.imagecrop[*] eq max(data.imagecrop) )
;         y0 = indmax[0] / imsize[0]
;         x0 = indmax[0] mod imsize[0]


        
;         parinfo=replicate({value:0.D, fixed:0, $
;                            limited:[0,0], limits:[0.D,0]}, 8)
        
       ;define parinfo, which gives limits to
       ;parameters as well as starting values
;         parinfo[*].value = [1.0, (gals.reff)[b,i], $
;                             (gals.sersic_n)[b,i], $
;                             (gals.ba_sdss)[i], $
;                             gals.ba_sdss[i] lt 0.5 ? -0.9 : 0.9, $
;                             x0, y0, 0.78 ]
        
;         parinfo[0].limited = [1,0]
;         parinfo[0].limits = [1D-8, 0]
;         parinfo[2].limited = [1,1]
;         parinfo[2].limits = [1D-8, 14.0] ; blanton uses 6.0
;         parinfo[1].limited = [1,0]
;         parinfo[1].limits = [1D-8,0]
;         parinfo[3].limited = [1,0]
;         parinfo[3].limits = [0.0D,1.0000001D]
;         parinfo[4].limited = [1,1]
;         parinfo[4].limits = [-2.D, 2.D]
;         parinfo[5].limited = [1,1]
;         parinfo[5].limits = [0.0, imsize[0]+0.05]
;         parinfo[6].limited = [1,1]
;         parinfo[6].limits = [0.0, imsize[1]+0.05]
        
;         start_params = parinfo[*].value
        

        ;normalize maximum pixel value to maximum pixel flux
;         norm = max(data.imagecrop) / $
;         max(pixelfluxpsf(x, y, start_params, _EXTRA={psfImage:data.psfcrop}))
        
;         parinfo[0].value = norm
;         start_params[0] = norm
                                ;DON'T USE errors, from DATA.ERR,
                                ;these are inverse variances, not
                                ;errors, rely on weights keyword for
                                ;errors!!!


;         temp = start_params


;         for fit=0, nfits do begin
        
;             again=1
;             fitstat=1
;             times = 0
            
;             while( again ne 0 and fitstat ne 5 and times lt 10) do begin
               
;                 start_params = temp
;                 again = 0
              
;                 if( fit ne 0 ) then begin
;                     params=mpfit2dfun('pixelfluxpsf', x, y, data.imagecrop, $
;                                       data.ivarcrop, $
;                                       start_params, parinfo=parinfo, $ 
;                                       weights=weights(data.ivarcrop, $
;                                                       data.imagecrop, $
;                                                       start_params, $
;                                                       gals.rpetro[i]*$
;                                                       fit_size[fit-1] ), $
;                                       status=fitstat, /quiet, $
;                                       perror=errors, covar=covar, $
;                                       dof=dof, bestnorm=bn,$
;                                       maxiter=100, $
;                                       functargs={psfImage:data.psfcrop})
;                 endif else begin
;                     params=mpfit2dfun('pixelfluxpsf', x, y, data.imagecrop, $
;                                       data.ivarcrop, $
;                                       start_params, parinfo=parinfo, $ 
;                                       weights=weights(data.ivarcrop, $
;                                                       data.imagecrop, $
;                                                       start_params, $
;                                                       0.0D ), $
;                                       status=fitstat, /quiet, $
;                                       perror=errors, covar=covar, $
;                                       dof=dof, bestnorm=bn,$
;                                       maxiter=100, $
;                                       functargs={psfImage:data.psfcrop})
;                 endelse
                
;                 if( params[3] gt 1.0 ) then begin
;                     again = 1
;                     params[3] = 1.0/params[3]
;                     params[7] -= !PI/2.0D
;                 endif
                
;                 temp = params
;                 times += 1
;             endwhile
            
;             if(fitstat eq 5 or times gt 9) then params[*] = 0.0
           
        ;write parameters, errors, and covariance to output array
        ; for the TOTAL fit only
 ;            if( fit eq 0 ) then begin
;                 output[i].DOF[b] = dof
;                 output[i].CHI_SQUARE[b] = bn
;                 output[i].SERSIC_NORM[b] = params[0]
;                 output[i].SERSIC_NORM_ERR[b] = errors[0]
;                 output[i].SERSIC_N[b] = params[2]
;                 output[i].SERSIC_N_ERR[b] = errors[2]
;                 output[i].SERSIC_R50[b] = params[1] * 0.396 ;to get arcsecs
;                 output[i].SERSIC_R50_ERR[b] = errors[1] * 0.396
;                 output[i].B_A[b] = params[3]
;                 output[i].B_A_ERR[b] = errors[3]
;                 output[i].C_ELLIP[b] = params[4]
;                 output[i].C_ELLIP_ERR[b] = errors[4]
;                 output[i].X_CENTER[b] = params[5] + xcrop
;                 ;in pixels with center indexed by integer
;                 output[i].X_CENTER_ERR[b] = errors[5]
;                 output[i].Y_CENTER[b] = params[6] + ycrop
;                 ;recenter to large image
;                 output[i].Y_CENTER_ERR[b] = errors[6]
;                                         ;phi in radians
;                 output[i].PHI_MAJ[b] = params[7]
;                 output[i].PHI_MAJ_ERR[b] = errors[7]
;                 output[i].COVAR[*,b] = covar
;                 output[i].SERSIC_FLUX[b] = totalsersicflux(params)
;             endif else begin
;                 output[i].FITSIZE[fit-1] = fit_size[fit-1]
;                 output[i].FS_SERSIC_N[fit-1,b] = params[2]
;                 output[i].FS_SERSIC_R50[fit-1,b] = params[1]*0.396
;                 output[i].FS_B_A[fit-1,b] = params[3]*0.396
;             endelse
;        endfor     
        
        ;output arrays for 

        
        ;display residuals
;          c = dindgen((size( data.image, /dimensions ))[0])
;          r = dindgen((size( data.image, /dimensions ))[1])
;         xdisp = c#(r*0+1)
;          ydisp = (c*0+1)#r
;         params[5] += xcrop
;          params[6] += ycrop
;        atv, data.imagecrop - pixelfluxpsf( x, y, params )
;        print, params

     
    endfor
endfor

;write output fits file: name of file=sersicHHhdecmindecmax.fits
rahrmin = floor( ramin / 360. * 24.0 )
outname = string( rahrmin, truncd(decmin,-1), roundd(decmax,-1), $
                  format='("/u/clackner/work/bulges/bulgeBranch/output/bd",i02,f+05.1,f+05.1,".fits")')
mwrfits, output, outname


END
