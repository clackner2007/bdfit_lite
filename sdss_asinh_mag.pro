;-------------------------------
; Claire Lackner
; sdss_asinh_mag
;
;  returns the asinh magnitude given the flux in maggies and the band

; output the magnitudes to the table
; these are SDSS asinh mags, as per the DR7 website
; m=-2.5/ln(10)*(asinh((f/f_0)/2b)+ln(b))
; where b is given for each band as:
;   b_asinh=[1.4e-10,0.9e-10,1.2e-10,1.8e-10,7.4e-10]

;---------------

FUNCTION sdss_asinh_mag, flux, filter=filter

b_asinh=[1.4e-10,0.9e-10,1.2e-10,1.8e-10,7.4e-10]
b=b_asinh[filter]

if(n_elements(b) gt 1) and n_elements(flux) ne n_elements(b) then begin
    fx = (size(flux, /dimensions))[0]
    fy = (size(flux, /dimensions))[1]
    b_new = dblarr(fx,fy)
    for i=0L, fx-1 do b_new[i,*] = b[i]
    b=b_new
    
endif

return, -2.5/alog(10.0)*(asinh(flux/(2.0*b))+alog(b))

END
