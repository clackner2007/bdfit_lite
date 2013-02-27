;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
; Jan 11 2010
; 
;  gets sample of galaxies given run, camcol, fiel,d id, rerun info
; put the sample into and output file called outfile, this is the prep
; for the fitter, as it contains most (NOT sersic info, or spectro
; info) of the info the fitter wants
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO getdata, galFile, outfile, vdisp=vdisp, z=z, min=min, max=max

;path='/peyton/scr/depot1/clackner/vagc2/'
;filephoto=path+'object_sdss_imaging.fits'
;filespectro=path+'object_sdss_spectro.fits'
;filesersic='/peyton/scr/depot1/clackner/vagc2/sersic/sersic_catalog.fits'


;photo=mrdfits( filephoto,1,rows=numrange[0] + indgen(numrange[1]-numrange[0],/long) )
;spectro=mrdfits(filespectro,1,rows=numrange[0] + indgen(numrange[1]-numrange[0],/long))
;sersic=mrdfits(filesersic,1,rows=numrange[0] + indgen(numrange[1]-numrange[0],/long))

galaxyList = mrdfits(galFile,1,range=[min,max])

properties = sdss_readobjlist(galaxyList.RUN, $
                              galaxyList.CAMCOL, $
                              galaxyList.FIELD, $
                              galaxyList.ID, $
                              rerun=galaxyList.RERUN, /silent,$
                              EXTRA={select_tags:['PIXSCALE','PSF_FWHM','AB_*', 'FRACPSF', $
                                                  '*FLUX*', '*_LNL', $
                                                  'R_*', 'ROWC', 'COLC', $
                                                  'RA', 'DEC', $
                                                  '*PETRO*', '*PHI*',$
                                                  'NMGYPERCOUNT', $
                                                  'EXTINCTION', 'NMGYPERCOUNT',$
                                                  'RESOLVE*', '*FLAG*', $
                                                  'NCHILD']})

;check the photometry flags and only keep the good ones
; see http://www.astro.princeton.edu/~rhl/flags.html and
; http://www.sdss.org/dr7/products/catalogs/flags.html for information
; on flags and selecting good phtometry

;NOTE: the FLAGS[0:4] are for the ugriz, colors, 2=g
;select for things that are not duplicates and not on the edge
i = (properties[*].RESOLVE_STATUS and $
      sdss_flagval('RESOLVE_STATUS','RUN_PRIMARY')) ne 0 $
  and (properties[*].RESOLVE_STATUS and $
       sdss_flagval('RESOLVE_STATUS','RUN_EDGE')) eq 0

;select for non-blended and not bright objects
j = ((properties[*].NCHILD eq 0 ) or $  ;no children
     (properties[*].FLAGS[2] and 2L^3) eq 0) and $ ;or not blended
  (properties[*].FLAGS[2] and 2L^6) eq 0 and $ ;and not deblended
  (properties[*].FLAGS[2] and 2L^1) eq 0 ;and not bright

;check there are no blending problems
DEBLEND_PROB1 = (2L^5 or 2L^19 ) ;PEAKCENTER or NOTCHECKED
DEBLEND_PROB2 = 2L^14 ;DEBLEND_NOPEAK
;interp problems, flag 2
INTERP_PROB2 = ( 2L^8 or 2L^15) ;PSF_FLUX_INTERP or BAD_COUNTS_ERROR
INTERP_CENTER2 = 2L^12
CR1 = 2L^12

k = (properties[*].FLAGS[2] and 2L^28) ne 0 and $ ;BINNED1
  (properties[*].FLAGS[2] and 2L^01) eq 0 and $ ;!BRIGHT
  (properties[*].FLAGS[2] and 2L^18) eq 0 and $ ;!SATUR(ated)
;     (properties[*].FLAGS[2] and 2L^02) eq 0 and $ ;!EDGE unset for now
((properties[*].FLAGS[2] and 2L^03) eq 0 or $ ; !BLENDED
 (properties[*].FLAGS[2] and 2L^06) ne 0) and $ ;or DEBLENDED
  (properties[*].FLAGS[2] and 2L^07) eq 0 and $ ;!NOPROFILE
  not(((properties[*].FLAGS2[2] and INTERP_PROB2) ne 0 or $
       ((properties[*].FLAGS2[2] and INTERP_CENTER2) ne 0 and $
        (properties[*].FLAGS[2] and CR1) ne 0))) and $
  not(((properties[*].FLAGS[2] and DEBLEND_PROB1) ne 0 or $
       ((properties[*].FLAGS2[2] and DEBLEND_PROB2) ne 0 ))) and $
                                        ;properties[*].             ;still need to make sure psf error is small
((properties[*].FLAGS[2] and 2L^0) eq 0) ;!DEBLENDED_AS_MOVING

shortlist = where(i and j and k, Nobjs, complement=badphoto)
print, 'removed ',n_elements(galaxyList)-Nobjs,' from set due to bad photometry'

mwrfits, galaxyList[badphoto], outfile+'_badphoto.fits'
galaxyList = galaxyList[shortlist]
properties = properties[shortList]



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
              EXTINCTION:dblarr(5), NMGYPERCOUNT:dblarr(5), $
              PSFFLUX:dblarr(5), PETROFLUX_IVAR:dblarr(5),$
              MODELFLUX_IVAR:dblarr(5), PSFFLUX_IVAR:dblarr(5) }

list = replicate( list_entry, Nobjs )
list[*].RA = properties.RA
list[*].DEC = properties.DEC
for i=0L, n_elements(shortlist)-1 do begin
    list[i].GALID = sdss_galid(list[i].RA,list[i].DEC)
endfor
list[*].RUN = galaxyList.RUN
list[*].RERUN = galaxyList.RERUN
list[*].CAMCOL = galaxyList.CAMCOL
list[*].FIELD = galaxyList.FIELD
list[*].ID = galaxyList.ID
list[*].BLANTON_N = 4.0                 ;default for now
list[*].BLANTON_FLUX = 0.0 ;default
list[*].BLANTON_REFF = 3.0 ; default
list[*].PSF_FWHM = properties.PSF_FWHM / properties.PIXSCALE
list[*].AB_DEV = properties.AB_DEV
list[*].AB_EXP = properties.AB_EXP
list[*].PETROFLUX = properties.PETROFLUX
list[*].MODELFLUX = properties.MODELFLUX ;flux for model fixed to r band
                                       ;exp or deV (depending on which
                                       ;fits better)
list[*].PETROFLUX_IVAR = properties.PETROFLUX_IVAR
list[*].MODELFLUX_IVAR = properties.MODELFLUX_IVAR
list[*].PSFFLUX = properties.PSFFLUX
list[*].PSFFLUX_IVAR = properties.PSFFLUX_IVAR
list[*].FRACDEV = properties.FRACPSF    ;fraction of flux from dVc
                                      ;profile in dVc+Exp profile
list[*].DEVFLUX = properties.DEVFLUX
list[*].EXPFLUX = properties.EXPFLUX
list[*].CMODELFLUX = properties.FRACPSF*properties.DEVFLUX + $
  (1.-properties.FRACPSF)*properties.EXPFLUX ;composite model flux
list[*].R_DEV = properties.R_DEV
list[*].R_EXP = properties.R_EXP
list[*].EXP_LNL = properties.EXP_LNL ;likelihood of exponential model
list[*].DEV_LNL = properties.DEV_LNL
list[*].PETROR50 = properties.PETROR50
list[*].PETROR90 = properties.PETROR90
list[*].PIXSCALE = properties.PIXSCALE
if keyword_set(vdisp) then begin
    list[*].VDISP = galaxyList.vdisp[0]
    list[*].VDISP_ERR = galaxyList.vdisp[1]
endif
if keyword_set(z) then list[*].Z = galaxyList.z
list[*].ROWC = properties.ROWC
list[*].COLC = properties.COLC
;change calibrated position angle to uncalibrated one
list[*].PHI_EXP = (properties.PHI_EXP_DEG + $
                   properties.PHI_OFFSET) * (!pi/180.0)
list[*].PHI_OFFSET_RAD = properties.PHI_OFFSET*!pi/180.0
list[*].EXTINCTION = properties.EXTINCTION
list[*].NMGYPERCOUNT = properties.NMGYPERCOUNT


mwrfits, list, outfile

end
