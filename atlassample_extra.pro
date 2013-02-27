;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
; Jan 11 2010
; 
;  retrieves a set of galaxies from the VAGC witha given redshift
;  and/or other constraints, which can then by extracted from atlas images

;  add terms needed to run calibrated objects through kcorrect code
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; returns the sing (-1 or +1) of a number, 0 is plus 1
; FUNCTION sgn, num
; neg = ( num lt 0 ) * (-1)
; pos = (num gt 0 or num eq 0)
; return, neg+pos
; end

;reads vagc catalog, to extract 
function atlassample_extra, numrange, zrange, rerun=rerun


path='/peyton/scr/depot1/clackner/vagc2/'
;path='/u/dss/redux/vagc/vagc0/sdss/'
filephoto=path+'object_sdss_imaging.fits'
filespectro=path+'object_sdss_spectro.fits'
filesersic='/peyton/scr/depot1/clackner/vagc2/sersic/sersic_catalog.fits'


photo=mrdfits( filephoto,1,rows=numrange[0] +indgen(numrange[1]-numrange[0],/long),$
             columns=['RA', 'DEC', 'RUN', 'RERUN', 'CAMCOL', 'FIELD',$
                      'ID', 'VAGC_SELECT', 'FLAGS', 'FLAGS2',$
                      'RESOLVE_STATUS', 'NCHILD', 'OBJC_FLAGS',$
                      'OBJC_FLAGS2','MODELFLUX'])
spectro=mrdfits(filespectro,1,rows=numrange[0] +indgen(numrange[1]-numrange[0],/long),$
                columns=['VDISP', 'VDISP_ERR', 'Z', 'SDSS_SPECTRO_TAG'])
sersic=mrdfits(filesersic,1,rows=numrange[0] +indgen(numrange[1]-numrange[0],/long),$
               columns=['SERSIC_N', 'SERSIC_FLUX', 'SERSIC_R50'])

;pass the main sample criterion
vagc = (photo[*].VAGC_SELECT and 4) ne 0

redshift = spectro[*].SDSS_SPECTRO_TAG ne -1 and $
  spectro[*].Z lt zrange[1] and  $
  spectro[*].Z gt zrange[0]

vagclist = where(vagc and redshift, nv)
if( nv eq 0) then begin
    print, "didn't find any in redshift range"
    return, -1
endif


;change the rerun number of list if desired
if keyword_set(rerun) then begin
    change_rerun, photo[vagclist].RUN, photo[vagclist].CAMCOL, photo[vagclist].FIELD, $
      photo[vagclist].RA, photo[vagclist].DEC, $
      photo[vagclist].MODELFLUX[2], $
      rerun=rerun, objid=new_id, properties=properties, $
      tags=['PIXSCALE','PSF_FWHM','AB_*', $
            'FRACPSF', $
            '*FLUX*', '*_LNL', $
            'R_*', 'ROWC', 'COLC',$
            '*PHI*','EXTINCTION', 'NMGYPERCOUNT',$
            'RESOLVE*','*FLAG*', 'NCHILD',$
            '*PETRO*']

    keep = where(new_id ne -1, nk)
    print, nk,'/ ',n_elements(vagclist)
    if( nk eq 0 ) then begin
        print, 'no new rerun objc found'
        return, -1
    endif
    
    vagclist = vagclist[where(new_id ne -1)]
    rerun = intarr(n_elements(vagclist)) + rerun
    properties = properties[where(new_id ne -1)]
    new_id = new_id[where(new_id ne -1)]
   
    
endif else begin
    new_id = photo[vagclist].ID
    rerun=photo[vagclist].RERUN

    properties = sdss_readobjlist(photo[vagclist].RUN, $
                              photo[vagclist].CAMCOL, $
                              photo[vagclist].FIELD, $
                              new_id, $
                              rerun=rerun, /silent,$
                              EXTRA={select_tags:['PIXSCALE','PSF_FWHM','AB_*', $
                                                  'FRACPSF', $
                                                  '*FLUX*', '*_LNL', $
                                                  'R_*', 'ROWC', 'COLC',$
                                                  '*PHI*','EXTINCTION', 'NMGYPERCOUNT',$
                                                  'RESOLVE*','*FLAG*', 'NCHILD',$
                                                  '*PETRO*']})
endelse

    

if n_elements(new_id) eq 0 then begin
    print, 'no vagc objects found'
    return, -1
endif

if not is_struct(properties) then begin
    print, "didn't find files"
    return, -1
endif

; see http://www.astro.princeton.edu/~rhl/flags.html and
; http://www.sdss.org/dr7/products/catalogs/flags.html for information
; on flags and selecting good phtometry

;NOTE: the FLAGS[0:4] are for the ugriz, colors, 2=r
;select for things that are not duplicates and not on the edge
i = (properties.RESOLVE_STATUS and $
      sdss_flagval('RESOLVE_STATUS','RUN_PRIMARY')) ne 0 $
  and (properties.RESOLVE_STATUS and $
       sdss_flagval('RESOLVE_STATUS','RUN_EDGE')) eq 0

;select for non-blended and not bright objects
j = (properties.NCHILD eq 0 ) and $     ;no children
  ((properties.OBJC_FLAGS and 2L^3) eq 0 or $ ;and not blended
   (properties.OBJC_FLAGS and 2L^6) ne 0) and $ ;or nodeblended
  (properties.OBJC_FLAGS and 2L^1) eq 0 ;and not bright

;check there are no blending problems
DEBLEND_PROB1 = (2L^5 or 2L^19 ) ;PEAKCENTER or NOTCHECKED
DEBLEND_PROB2 = 2L^14 ;DEBLEND_NOPEAK
;interp problems, flag 2
INTERP_PROB2 = ( 2L^8 or 2L^15) ;PSF_FLUX_INTERP or BAD_COUNTS_ERROR
INTERP_CENTER2 = 2L^12
CR1 = 2L^12

k = (((properties.OBJC_FLAGS and 2L^28) ne 0) and $ ;BINNED1
     ((properties.FLAGS[2] and 2L^18) eq 0) and $ ;!SATUR(ated)
;     (properties.FLAGS[2] and 2L^02) eq 0 and $ ;!EDGE unset for now
     ((properties.FLAGS[2] and 2L^07) eq 0) and $ ;!NOPROFILE
       not((((properties.FLAGS2[2] and INTERP_PROB2) ne 0) or $
            (((properties.FLAGS2[2] and INTERP_CENTER2) ne 0) and $
             ((properties.FLAGS[2] and CR1) ne 0)))) and $
       not((((properties.FLAGS[2] and DEBLEND_PROB1) ne 0)  $
            or (((properties.FLAGS2[2] and DEBLEND_PROB2) ne 0)))) $
       and $
       ((properties.FLAGS[2] and 2L^0) eq 0)) ;!DEBLENDED_AS_MOVING

shortlist = where(i and j and k, Nobjs )
photo_shortlist = vagclist[shortlist]

if( Nobjs eq 0 ) then begin
    print, 'no objects found'
    return, -1
endif

print, 'original numbers '
print, n_elements(vagclist)
print, 'new numbers '
print, Nobjs, n_elements(photo_shortlist)

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

list = replicate( list_entry, n_elements(shortlist) )

list[*].RA = photo[photo_shortlist].RA
list[*].DEC = photo[photo_shortlist].DEC
for i=0L, n_elements(shortlist)-1 do begin
    list[i].GALID = sdss_galid(list[i].RA,list[i].DEC)
    endfor
list[*].RUN = photo[photo_shortlist].RUN
list[*].RERUN = rerun[shortlist]
list[*].CAMCOL = photo[photo_shortlist].CAMCOL
list[*].FIELD = photo[photo_shortlist].FIELD
list[*].ID = new_id[shortlist]
list[*].BLANTON_N = sersic[photo_shortlist].SERSIC_N
list[*].BLANTON_FLUX = sersic[photo_shortlist].SERSIC_FLUX
list[*].BLANTON_REFF = sersic[photo_shortlist].SERSIC_R50 / $
  properties[shortlist].PIXSCALE        ;photo[shortlist].PIXSCALE
list[*].PSF_FWHM = properties[shortlist].PSF_FWHM / properties[shortlist].PIXSCALE
list[*].AB_DEV = properties[shortlist].AB_DEV
list[*].AB_EXP = properties[shortlist].AB_EXP
list[*].PETROFLUX = properties[shortlist].PETROFLUX
list[*].MODELFLUX = properties[shortlist].MODELFLUX ;flux for model fixed to r band
                                       ;exp or deV (depending on which
                                       ;fits better)
list[*].PETROFLUX_IVAR = properties[shortlist].PETROFLUX_IVAR
list[*].MODELFLUX_IVAR = properties[shortlist].MODELFLUX_IVAR
list[*].PSFFLUX = properties[shortlist].PSFFLUX
list[*].PSFFLUX_IVAR = properties[shortlist].PSFFLUX_IVAR

list[*].FRACDEV = properties[shortlist].FRACPSF      ;fraction of flux from dVc
                                      ;profile in dVc+Exp profile
list[*].DEVFLUX = properties[shortlist].DEVFLUX
list[*].EXPFLUX = properties[shortlist].EXPFLUX
list[*].CMODELFLUX = properties[shortlist].FRACPSF*properties[shortlist].DEVFLUX + $
  (1.-properties[shortlist].FRACPSF)*properties[shortlist].EXPFLUX ;composite model flux
list[*].R_DEV = properties[shortlist].R_DEV
list[*].R_EXP = properties[shortlist].R_EXP
list[*].EXP_LNL = properties[shortlist].EXP_LNL ;likelihood of exponential model
list[*].DEV_LNL = properties[shortlist].DEV_LNL
list[*].PETROR50 = properties[shortlist].PETROR50
list[*].PETROR90 = properties[shortlist].PETROR90
list[*].PIXSCALE = properties[shortlist].PIXSCALE
list[*].VDISP = spectro[photo_shortlist].VDISP
list[*].VDISP_ERR = spectro[photo_shortlist].VDISP_ERR
list[*].Z = spectro[photo_shortlist].Z
list[*].ROWC = properties[shortlist].ROWC
list[*].COLC = properties[shortlist].COLC
;put into uncalibrated coordinates
list[*].PHI_EXP = (properties[shortlist].PHI_EXP_DEG +$
                   properties[shortlist].PHI_OFFSET) * (!pi/180.0)
list[*].PHI_OFFSET_RAD = properties[shortlist].PHI_OFFSET*!pi/180.0
list[*].EXTINCTION = properties[shortlist].EXTINCTION
list[*].NMGYPERCOUNT = properties[shortlist].NMGYPERCOUNT


;check image, comment out for actual code
if 0 then begin
    for jj=0L, n_elements(shortlist)-1 do begin
        s1=sdss_atlas_image(photo[photo_shortlist[jj]].RUN, photo[photo_shortlist[jj]].CAMCOL,$
                            photo[photo_shortlist[jj]].FIELD,photo[photo_shortlist[jj]].ID,$
                            rerun=137, filter=2)
        s2=sdss_atlas_image(photo[photo_shortlist[jj]].RUN, photo[photo_shortlist[jj]].CAMCOL,$
                            photo[photo_shortlist[jj]].FIELD,new_id[shortlist[jj]],$
                            rerun=301, filter=2)
        if n_elements(s1) le 1 or n_elements(s2) le 1 then continue
        s=dblarr(max([(size(s1,/dimensions))[0],(size(s2,/dimensions))[0]]),$
                 (size(s1,/dimensions))[1]+(size(s2,/dimensions))[1])
        s[0,0]=s2-1000.
        s[0,(size(s2,/dimensions))[1]]=s1-1000.
        print, sdss_galid(photo[photo_shortlist[jj]].RA,$
                          photo[photo_shortlist[jj]].DEC)
        atv,s
    endfor
endif

return, list

end
