;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
; Jan 11 2010
; 
;  retrieves a set of galaxies from the VAGC witha given redshift
;  and/or other constraints, which can then by extracted from atlas images
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; returns the sing (-1 or +1) of a number, 0 is plus 1
; FUNCTION sgn, num
; neg = ( num lt 0 ) * (-1)
; pos = (num gt 0 or num eq 0)
; return, neg+pos
; end

;reads vagc catalog, to extract 
function atlassample, numrange, zrange


path='/peyton/scr/depot1/clackner/vagc2/'
;path='/u/dss/redux/vagc/vagc0/sdss/'
filephoto=path+'object_sdss_imaging.fits'
filespectro=path+'object_sdss_spectro.fits'
filesersic='/peyton/scr/depot1/clackner/vagc2/sersic/sersic_catalog.fits'


photo=mrdfits( filephoto,1,rows=numrange[0] + indgen(numrange[1]-numrange[0],/long) )
spectro=mrdfits(filespectro,1,rows=numrange[0] + indgen(numrange[1]-numrange[0],/long))
sersic=mrdfits(filesersic,1,rows=numrange[0] + indgen(numrange[1]-numrange[0],/long))

;pass the main sample criterion
vagc = (photo[*].VAGC_SELECT and 4) ne 0

; see http://www.astro.princeton.edu/~rhl/flags.html and
; http://www.sdss.org/dr7/products/catalogs/flags.html for information
; on flags and selecting good phtometry

;NOTE: the FLAGS[0:4] are for the ugriz, colors, 2=g
;select for things that are not duplicates and not on the edge
i = (photo[*].RESOLVE_STATUS and $
      sdss_flagval('RESOLVE_STATUS','RUN_PRIMARY')) ne 0 $
   and (photo[*].RESOLVE_STATUS and $
        sdss_flagval('RESOLVE_STATUS','RUN_EDGE')) eq 0

;select for non-blended and not bright objects
j = ((photo[*].NCHILD eq 0 ) or $ ;no children
   (photo[*].FLAGS[2] and 2L^3) eq 0) and $ ;or not blended
   (photo[*].FLAGS[2] and 2L^6) eq 0 and $ ;and not deblended
   (photo[*].FLAGS[2] and 2L^1) eq 0   ;and not bright

;check there are no blending problems
DEBLEND_PROB1 = (2L^5 or 2L^19 ) ;PEAKCENTER or NOTCHECKED
DEBLEND_PROB2 = 2L^14 ;DEBLEND_NOPEAK
;interp problems, flag 2
INTERP_PROB2 = ( 2L^8 or 2L^15) ;PSF_FLUX_INTERP or BAD_COUNTS_ERROR
INTERP_CENTER2 = 2L^12
CR1 = 2L^12

k = ((photo[*].FLAGS[2] and 2L^28) ne 0 and $ ;BINNED1
     (photo[*].FLAGS[2] and 2L^01) eq 0 and $ ;!BRIGHT
     (photo[*].FLAGS[2] and 2L^18) eq 0 and $ ;!SATUR(ated)
;     (photo[*].FLAGS[2] and 2L^02) eq 0 and $ ;!EDGE unset for now
     ((photo[*].FLAGS[2] and 2L^03) eq 0 or $ ; !BLENDED
      (photo[*].FLAGS[2] and 2L^06) ne 0) and $ ;or DEBLENDED
     (photo[*].FLAGS[2] and 2L^07) eq 0 and $ ;!NOPROFILE
     not(((photo[*].FLAGS2[2] and INTERP_PROB2) ne 0 or $
          ((photo[*].FLAGS2[2] and INTERP_CENTER2) ne 0 and $
           (photo[*].FLAGS[2] and CR1) ne 0))) and $
     not(((photo[*].FLAGS[2] and DEBLEND_PROB1) ne 0 or $
          ((photo[*].FLAGS2[2] and DEBLEND_PROB2) ne 0 ))) and $
                                        ;photo[*].             ;still need to make sure psf error is small
      (photo[*].FLAGS[2] and 2L^0) eq 0 )     ;!DEBLENDED_AS_MOVING

redshift = spectro[*].SDSS_SPECTRO_TAG ne -1 and $
  spectro[*].Z lt zrange[1] and  $
  spectro[*].Z gt zrange[0]

shortlist = where( vagc and i and j and k and redshift, Nobjs )

if( Nobjs eq 0 ) then begin
    print, 'no objects found'
    return, -1
endif

moredata = sdss_readobjlist(photo[shortlist].RUN, $
                            photo[shortlist].CAMCOL, $
                            photo[shortlist].FIELD, $
                            photo[shortlist].ID, $
                            rerun=photo[shortlist].RERUN, /silent,$
                            EXTRA={select_tags:['PIXSCALE','PSF_FWHM','AB_*', $
                                                'FRACPSF', $
                                                '*FLUX*', '*_LNL', $
                                                'R_*', 'ROWC', 'COLC',$
                                                '*PHI*','EXTINCTION', 'NMGYPERCOUNT']})

list_entry = {GALID:'SDSSJhhmmss.ss+ddmmss.s',$
              RA:0.0D, DEC:0.0D, RUN:1L, CAMCOL:1L, FIELD:1L, ID:1L, $
              RERUN:1L, BLANTON_N:dblarr(5), BLANTON_FLUX:dblarr(5), $
              BLANTON_REFF:dblarr(5), PSF_FWHM:dblarr(5), $
              AB_DEV:dblarr(5), AB_EXP:dblarr(5), $
              PETROFLUX:dblarr(5), MODELFLUX:dblarr(5), $
              FRACDEV:dblarr(5), DEVFLUX:dblarr(5), $
              EXPFLUX:dblarr(5), CMODELFLUX:dblarr(5), $
              R_DEV:dblarr(5), R_EXP:dblarr(5), $
              EXP_LNL:dblarr(5), DEV_LNL:dblarr(5), $
              PETROR50:dblarr(5), PETROR90:dblarr(5), $
              PIXSCALE:dblarr(5), $
              VDISP:0.0D,VDISP_ERR:0.0D, Z:0.0D, $
              ROWC:dblarr(5), COLC:dblarr(5), $
              PHI_EXP:dblarr(5), PHI_OFFSET_RAD:dblarr(5),$
              EXTINCTION:dblarr(5), NMGYPERCOUNT:dblarr(5) }

list = replicate( list_entry, n_elements(shortlist) )

list[*].RA = photo[shortlist].RA
list[*].DEC = photo[shortlist].DEC
for i=0L, n_elements(shortlist)-1 do begin
    list[i].GALID = sdss_galid(list[i].RA,list[i].DEC)
    endfor
list[*].RUN = photo[shortlist].RUN
list[*].RERUN = photo[shortlist].RERUN
list[*].CAMCOL = photo[shortlist].CAMCOL
list[*].FIELD = photo[shortlist].FIELD
list[*].ID = photo[shortlist].ID
list[*].BLANTON_N = sersic[shortlist].SERSIC_N
list[*].BLANTON_FLUX = sersic[shortlist].SERSIC_FLUX
list[*].BLANTON_REFF = sersic[shortlist].SERSIC_R50 / moredata.PIXSCALE ;photo[shortlist].PIXSCALE
list[*].PSF_FWHM = moredata.PSF_FWHM / moredata.PIXSCALE
list[*].AB_DEV = moredata.AB_DEV
list[*].AB_EXP = moredata.AB_EXP
list[*].PETROFLUX = moredata.PETROFLUX
list[*].MODELFLUX = moredata.MODELFLUX ;flux for model fixed to r band
                                       ;exp or deV (depending on which
                                       ;fits better)
list[*].FRACDEV = moredata.FRACPSF    ;fraction of flux from dVc
                                      ;profile in dVc+Exp profile
list[*].DEVFLUX = moredata.DEVFLUX
list[*].EXPFLUX = moredata.EXPFLUX
list[*].CMODELFLUX = moredata.FRACPSF*moredata.DEVFLUX + $
  (1.-moredata.FRACPSF)*moredata.EXPFLUX ;composite model flux
list[*].R_DEV = moredata.R_DEV
list[*].R_EXP = moredata.R_EXP
list[*].EXP_LNL = moredata.EXP_LNL ;likelihood of exponential model
list[*].DEV_LNL = moredata.DEV_LNL
list[*].PETROR50 = photo[shortlist].PETROR50
list[*].PETROR90 = photo[shortlist].PETROR90
list[*].PIXSCALE = moredata.PIXSCALE
list[*].VDISP = spectro[shortlist].VDISP
list[*].VDISP_ERR = spectro[shortlist].VDISP
list[*].Z = spectro[shortlist].Z
list[*].ROWC = moredata.ROWC
list[*].COLC = moredata.COLC
;put into uncalibrated coordinates
list[*].PHI_EXP = (moredata.PHI_EXP_DEG +$
                   moredata.PHI_OFFSET) * (!pi/180.0)
list[*].PHI_OFFSET_RAD = moredata.PHI_OFFSET*!pi/180.0
list[*].EXTINCTION = moredata.EXTINCTION
list[*].NMGYPERCOUNT = moredata.NMGYPERCOUNT

return, list

end
