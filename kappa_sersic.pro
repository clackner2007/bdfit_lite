;____________________________________
;
; Claire Lackner
;  Sept 7, 2009
;
;  kappa_sersic.pro
;  
; function for computing facotr of exponent in Sersic profile
; such that half the light is contained within Reff
;____________________________________

FUNCTION kappa_sersic, n
if( n_elements(n) eq 1 ) then begin
    if n lt 1.0e-18 then return, 0.0 $
      else return, n*exp(0.6950-0.1789/n)
endif else begin
    k = n
    k[where(n gt 1.0e-18)] *= exp(0.6950-0.1789/n[where(n gt 1.0e-18)])
;approximation from Lima Neta, 1999
;ee=0.6950-0.1789/n
;k = n*exp(ee)
    return, k
endelse

END
