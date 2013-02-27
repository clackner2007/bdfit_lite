;computes the bulge to total error given the covariance martix and the
;fit parameters

FUNCTION bulgetotot_error, params, covar, cutoff=cutoff

if( not keyword_set(cutoff)) then cutoff = 0

;parameters needed to calculate bulge to total
bttparams = [0,1,3,8,9,11]

bttcovar = (covar[[bttparams],*])[*,[bttparams]]

h = params[bttparams]*0.00005
zero=where(h eq 0.0,z)
if(z gt 0) then h[zero] = 0.00005
derivs = dblarr(6)
for i=0, 5 do begin
    delta = dblarr(n_elements(params))
    delta[bttparams[i]] = h[i]
    derivs[i] = (bulgetotot(params+delta, cutoff=cutoff) - $
                 bulgetotot(params, cutoff=cutoff)) / $
      delta[bttparams[i]]
endfor

sigma2_btt = total(derivs*(bttcovar#derivs))

;multiplicative factors of disk and bulge flux
; func = bulgetotot(params, cutoff=1)

; beta = totalsersicflux([params[6:10],params[12:13],params[11]])/$
;   (params[ib]*params[rb]^2*params[qb])
; delta = totalsersicflux([params[0:4],params[12:13],params[5]])/$
;   (params[id]*params[rd]^2*params[qd])

return, sqrt(sigma2_btt)

END

