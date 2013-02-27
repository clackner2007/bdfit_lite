;----------------------------------
; Claire Lackner
; Mar 10 2010
;
;
; bulgetotot.pro
;
; Calculates the bulge to toal flux ratio for a set of 14 bulge+disk
; parameters, currently assumes no cutoff which is nonsense, but
; that's what the current fits do as well
;
;
;---------------------------------------------

FUNCTION bulgetotot, params, cutoff=cutoff

pbulge = [params[8:15]]
pdisk = [params[0:7]]
if( pbulge[2] lt 1.0e-8 and pdisk[2] lt 1.0e-8 ) then return, -1.
if( not keyword_set(cutoff)) then cutoff = 0

fluxbulge = totalsersicflux( pbulge, cutoff=cutoff)
btt = fluxbulge/(fluxbulge+totalsersicflux(pdisk, cutoff=cutoff))
return, btt

END
