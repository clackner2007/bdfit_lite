;computers the error in the totalsersicflux

FUNCTION totalsersicflux_error, params, covar, disk=disk, bulge=bulge, cutoff=cutoff

if keyword_set(bulge) then begin
    if not keyword_set(disk) then begin
        vc_par = [8,9,11]
        v_par = [0,1,3]
        my_params = dblarr(8)
        my_params[0:4] = params[8:12]
        my_params[5:6] = params[13:14]
        my_params[7] = params[15]
        n=3
    endif else begin
        v_par = [0,1,3,8,9,11]
        vc_par = v_par
        n=6
        my_params=params
    endelse
endif else begin
    v_par = [0,1,2,3]
    vc_par = v_par
    if n_elements(params) gt 8 then begin
        my_params = dblarr(8)
        my_params[0:4] = params[0:4]
        my_params[5:6] = params[5:6]
        my_params[7] = params[7]
    endif else my_params=params
    n=4
endelse

my_covar = (covar[[vc_par],*])[*,[vc_par]]

dell=1.e-6
h = my_params[v_par]*dell
zero=where(h eq 0.0,z)
if(z gt 0) then h[zero] = dell
derivs = dblarr(n)
for i=0, n-1 do begin
    delta = dblarr(n_elements(my_params))
    delta[v_par[i]] = h[i]
    derivs[i] = (totalsersicflux(my_params+delta, cutoff=cutoff) - $
                 totalsersicflux(my_params, cutoff=cutoff)) / $
      delta[v_par[i]]
endfor

sigma2_tsf = total(derivs*(my_covar#derivs))


return, sqrt(sigma2_tsf)

END
