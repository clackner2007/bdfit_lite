;+
;default_init_cond
;
; function to return default initial conditions for a fixed number of
;profiles
;
;ARGUMENTS=nprof, image, sersics=sersics, rescale=rescale
;
;if sersics isn't given the default is [4.0,1.0,1.0,...]
;
;RETURNS: {params: [8*nprof, 1d array of initial conditions],
;          fixed: [8*nprof, 1d array of fixed params, only C is fixed now]}
;-


FUNCTION default_init_cond, nprof, image, sersics=sersics, rescale=rescale
  params = dblarr(8*nprof)
  fixed_params = intarr(8*nprof)
  plist = indgen(nprof)*8
  
  if n_elements(sersics) eq 0 then begin
     sersics = dblarr(nprof) + 1.0
     sersics[0] = 4.0 
  endif
  
  ;fix the C values

  fixed_params[plist+4] = 1
  xs = getxy_start(image)
  x0 = xs[0]
  y0 = xs[1]
  
  params[0:7] = [1.0,10.,sersics[0], 0.7, 0.0, x0, y0, 0.2]

  for ip=1, nprof-1 do begin
     if keyword_set(rescale) then rad = 1.9^ip else rad = 1.9^ip*params[1] 
     params[ip*8:ip*8+7] = [1.0, rad, sersics[ip], 0.85, $
                            0.0, x0, y0, 0.2]
  endfor
  
  setscaling, image, params, fixed_params, rescale=rescale, fracs=dblarr(nprof)+1.0/nprof

  return, {params:params, fixed:fixed_params}
END
