;____________________________________
;
; Claire Lackner
;  Sept 7, 2009
;
;  totalsersicflux.pro
;
;  IDL function that returns the total sersic flux for the parameters:
;  a_norm      normalization of surface brightness profile
;  r_eff       effective radius (in pixels)
;  n_sersic    sersic index
;  q_ba        ellipticity (minor axis/major axis)
;  c_ellip     diskiness/boxiness parameter c=0 for ellipse
;              c > 0 for box, c < 0 for disk
;  x_0         central position in x
;  y_0         central position in y
;  phi         rotation angle of major axis
;____________________________________


;function for computing facotr of exponent in Sersic profile
; such that half the light is contained within Reff
FUNCTION ksers, n
;approximation from Lima Neta, 1999
k = n*exp(0.6950-0.1789/n)
return, k
END

;calculates the total flux of a sersic profile using the parameters p
;given in the order listed above
FUNCTION totalsersicflux, param, cutoff=cutoff

if(n_elements(param) ge 16 ) then begin
    flux = totalsersicflux([param[0:4],0.0,0.0,0.0], cutoff=cutoff)+$
      totalsersicflux([param[8:12],0.0,0.0,0.0],cutoff=cutoff)
    return, flux
endif


k = ksers(param[2])
n = param[2]
reff = param[1]
c2 = param[4]+2.
q = param[3]
Rc = !pi*c2 / (4.*beta(1./c2, 1.+1./c2, /double))

flux = 2.*!pi*reff^2 * param[0] * exp(k) * $
  n * k^(-2.*n)*gamma(2.*n)*q / Rc

; if the cutoff is set, assume it is set to the SDSS values 7-8 for a
; n>1 and 3-4 for an n<=1 profile from curve fitting of total flux,
; the following polynomail fits work:
;  for n>1 (mailny n=4) case:
;      cutoff_flux = totalflux*( 0.936587 - 2.7233496e-2*(n-4.) +
;      1.4125945e-3*(n-4.)^2 + 4.06329317e-4*(n-4.)^3
;      -1.82907176e-4*(n-4.)^4 +2.2355412e-5*(n-4.)^5 )
;  for 0.6<n<=1 (mainly n=1) case:
;      cutoff_flux = totalflux*(0.98042454e-1 - 0.06495062*(n-1.) -
;      0.02627548*(n-1.)^2 + 0.02277057*(n-1.)^3 - 0.21949023*(n-1)^4
;      - 0.34741306*(n-1.)^5)

if( keyword_set(cutoff) ) then begin
    if( (n - 1.0) lt 1.0e-5 ) then begin
        if( n eq 1.0 ) then flux *= 0.980472 else $
          flux *= 0.98042454 - 0.06495062*(n-1.) - $
          0.02627548*(n-1.)^2 + 0.02277057*(n-1.)^3 - $
          0.21949023*(n-1)^4 - 0.34741306*(n-1.)^5        
    endif else begin
        if( n eq 4.0 ) then flux *= 0.93666 else $
          flux *= 0.936587 - 2.7233496e-2*(n-4.) + $
          1.4125945e-3*(n-4.)^2 + 4.06329317e-4*(n-4.)^3 - $
          1.82907176e-4*(n-4.)^4 +2.2355412e-5*(n-4.)^5
    endelse
endif

return, flux
END
