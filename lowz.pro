;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
; Sept 15, 2009
; code to return list of lowz objects that satisfy given requirements
; in position space, the 

;returns a structur containing, object name (JHHMMSS.ss+DDMMSS.s),
;object sersic index from sdss, object 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;returns the sing (-1 or +1) of a number, 0 is plus 1
FUNCTION sgn, num
neg = ( num lt 0 ) * (-1)
pos = (num gt 0 or num eq 0)
return, neg+pos
end

;reads the sdss dr4 lowz catalog, i nthe given ra and dec range, dec
;given in degrees, and ra given in degrees
function lowz, ra_range, dec_range, status=status

status=1
;change to dr6 to match up with the images
path='/peyton/scr/depot1/clackner/vagc2/lowz_dr6/'
file=path+'lowz_plus_ned.dr6.fits'
path4 ='/peyton/scr/depot1/clackner/vagc2/lowz/'
file4 = path4+'lowz_catalog.dr4.fits'

fileextra=path+'lowz_measure.dr6.fits'

lz=mrdfits(file,1, /silent)
extra=mrdfits(fileextra,1,/silent)

hrdeg = where( lz[*].RA lt ra_range[1] and lz[*].RA gt ra_range[0] and $
               lz[*].DEC lt dec_range[1] and lz[*].DEC gt dec_range[0] and $
               (lz[*].QUALITY eq 0 or lz[*].QUALITY eq 1 ) and $
               lz[*].CLASS eq 'GALAXY' and $
               (lz[*].VAGC_SELECT and 4) ne 0, $ 
               number)

;check to make sure something is there
if( number eq 0 ) then begin
    status = -1
    return, -1
endif

;the right ascension is always positve
hrRA = floor(lz[hrdeg].RA/360.*24.)
minRA = floor((hrRA - (lz[hrdeg].RA/360.*24.))*60.*(-1.))
secRA = truncd(((lz[hrdeg].RA/360.*24. - hrRA)*60. - minRA)*60., (-2))

;need to worry about signs
degDEC = floor( abs(lz[hrdeg].DEC) )*sgn( lz[hrdeg].DEC )
minDEC = floor(abs(lz[hrdeg].DEC - degDEC)*60.)
secDEC = truncd((abs(lz[hrdeg].DEC-degDEC)*60. - minDEC)*60., (-1))


name=strarr( n_elements(hrdeg) )
n = fltarr( n_elements(hrdeg) )
reff = fltarr( n_elements(hrdeg) )
fld = strarr( n_elements(hrdeg) )
;writes object name, sersic index, the sersic effective radius (in pixels)
for i=0,(n_elements(hrdeg)-1) do begin
    index = hrdeg[i]
    posneg = string( sgn(lz[index].DEC) gt 0? '+': '-' )
    name[i] = string(hrRA[i], minRA[i], secRA[i], posneg, abs(degDEC[i]), $
                     minDEC[i], secDEC[i], $
                     format='("J",i02,i02,f05.2,a,i02,i02,f04.1)' )
    fld[i] = string( hrRA[i], format='(i02,"h")' )

;    n[i] = (lz[index].SERSIC_N)[2]
;    reff[i] = (lz[index].SERSIC_TH50)[2]*2.525


endfor

list = { name:name, folder:fld, sersic_n:lz[hrdeg].SERSIC_N, $
         reff:lz[hrdeg].SERSIC_TH50 * 2.525, flux:lz[hrdeg].SERSIC_FLUX, $
         run:lz[hrdeg].RUN, rerun:lz[hrdeg].RERUN, $
         camcol:lz[hrdeg].CAMCOL, id:lz[hrdeg].ID, field:lz[hrdeg].FIELD,$
         ra:lz[hrdeg].RA, dec:lz[hrdeg].DEC, $
         psf_fwhm:lz[hrdeg].PSF_FWHM*2.525, rpetro:extra[hrdeg].PETRORAD, $
         ba_sdss:extra[hrdeg].BA90 }

status = 1

return, list

end
