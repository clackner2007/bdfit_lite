;+
;____________________________________
;
; Claire Lackner
;  Sept 7, 2009
;
;  sersic.pro
;
;  IDL function that returns a sersic profile
; with the 8 free parameters, givin in this order:
;  a_norm      normalization of surface brightness profile
;  r_eff       effective radius (in pixels)
;  n_sersic    sersic index
;  q_ba        ellipticity (minor axis/major axis)
;  c_ellip     diskiness/boxiness parameter c=0 for ellipse
;              c > 0 for box, c < 0 for disk
;  x_0         central position in x
;  y_0         central position in y
;  phi         rotation angle of major axis
;
;
;the flux at a point x, y is:
;
;I * exp {-k(n) * [(r / R)^(1/n) - 1]}
;r = sqrt{ [ (x-x0)*cos(phi) + (y-y0)*sin(phi) ]^2 + 
;    [ ( (y-y0)*cos(phi) - (x-x0)*sin(phi) )/q ]^2 }
;
;where
; k(n) = exp[0.6950-0.1789/n]
;I take the functional form for k(n) from Lima Neto, 1999; but any reasonable one should work.
;
;____________________________________
;-

;function for computing radius of generalized ellipse
; param contains necessary parameters, and xy is the point at which 
; r is computed
FUNCTION radius, param, x, y
iq = 1.0/param[0]
c = param[1]
x0 = param[2]
y0 = param[3]
phi = param[4]

axtrans = abs((x - x0) * cos(phi) + (y - y0) * sin(phi))
aytransq = abs(((y - y0) * cos(phi) - (x - x0) * sin(phi))*iq)
if(c > 1.0e-12 ) then begin
    expo = c+2.0
    expo_inv = 1.0/expo
    r = ( axtrans^(expo) + aytransq^(expo) ) ^ expo_inv
endif else r = sqrt(axtrans^2 + aytransq^2)

return, r 
END

;has c=0, so don't need to do float exponents
FUNCTION radiusfast, param, x, y
iq = 1.0/param[0]

xx0 = x-param[2]
yy0 = y-param[3]
cp = cos(param[4])
sp = sin(param[4])

axtrans = (xx0) * cp + (yy0) * sp
aytransq = ((yy0) * cp - (xx0) * sp)*iq
r = sqrt(axtrans^2 + aytransq^2)

return, r

END

;sersic function, x is the point at which it is computed
; p is the set of parameters, given in the order
; listed above
; the extra keyword determines whether or not to use the cutoff radii
FUNCTION sersic, x, y, p, cutoff=cutoff

; attempt to make initial p[1] circularized effective radius
; (not along major axis)
p[1] = p[1] * sqrt(1.0/p[3])

if size(x, /n_dimensions) gt 1 then begin
    density = dblarr( (size( x, /dimensions))[0], $
                      (size( y, /dimensions))[1] )
endif else begin
    density= dblarr(n_elements(x))
endelse

if( abs(p[0]) lt 1.0e-18 or p[1] lt 1.0e-12 ) then return, density

if( abs(p[4]) lt 1.0e-12 ) then r=radiusfast( p[3:*],x,y ) $
else r = radius( p[3:*], x, y )
n = p[2]
k = kappa_sersic( n )
ir0 = 1.0/p[1]

;check if n is an easy number to work with
if( abs(n-1.0) lt 1.0e-12 ) then exponent = (-k) * ( r*ir0 - 1.0 ) $
else if( abs(n-4.0) lt 1.0e-12 ) then $
  exponent = (-k)*(sqrt(sqrt(r*ir0))-1.0) $
else if( abs(n-2.0) lt 1.0e-12 ) then $
  exponent = (-k)*(sqrt(r*ir0)-1.0) $
else exponent = (-k) * ( (r*ir0)^(1./n) - 1.0 )

;only compute exponent above lower bound, otherwise, leave value as
;zero
notzero = where( exponent gt -30, count )
if( count gt 0 ) then density[notzero] = p[0]*exp(exponent[notzero])

;cut the profile off
; smoothly past a certain size
;use the SDSS defined cutoffs

if( keyword_set(cutoff) ) then begin
    if n le 1.05 then begin
        fade = 3.*p[1]
        cut = 4.*p[1]
    endif else begin
        fade = 7.*p[1]
        cut = 8.*p[1]
    endelse
    denom = 1.0/(cut-fade)^2
    
    if size(x,/n_elements) ge 1 then begin
        truncate = where( (r gt fade) and (r le cut), t )
        cutoff = where( r gt cut, c )
        if( t ne 0 ) then $
          density[truncate] *= (1.-((r[truncate]- fade)^2*denom))^2
        if( c ne 0 ) then $
          density[cutoff] = 0.0
    endif else begin
        if(r gt cut ) then density = 0.0
        if( r gt fade ) then density*=(1.0-(r-fade)^2*denom)^2
    endelse
    
endif
return, density
END
