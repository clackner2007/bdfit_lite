;+
;sets the scaling of the profiles for the initial conditions
;
;-
PRO setscaling, image, params, fixed_params, rescale=rescale, fracs=fracs
  
  plist = indgen(n_elements(params)/8)*8 
  rescale = keyword_set(rescale)
  tot = max([total(image),10])
  nprof = n_elements(plist)
  
  totp = dblarr(n_elements(plist))
  for ip=0, n_elements(plist)-1 do begin
     if ((ip gt 0) and (rescale ne 0)) then rad = params[plist[ip]+1]*params[1] $
        else $
           rad = params[plist[ip]+1]
     totp[ip] = max([totalsersicflux([params[plist[ip]], $
                                     rad, $
                                     params[plist[ip]+2:plist[ip]+7]]), 1.0e-3])
  endfor

  if n_elements(fracs) eq 0 then fracs = dblarr(nprof)+1.0/nprof
  for ip=0, n_elements(plist)-1 do begin
     if fixed_params[plist[ip]] ne 0 then continue
     params[plist[ip]] *= fracs[ip]*tot/totp[ip]
  endfor
  
END
