;+
;NAME: init_cond
;
;AUTHOR:
; Claire Lackner
;
;DATE:
;  Jan 21 2014
;
;PURPOSE:
; Generate initial conditions for fitter, assuming some set of fixed
; values
;
;CALL:
; init_cond, params, fixed_params, image
; params = params, array with the number of params (8x number of profiles)
; fixed_params = mask with fixed parameters set to 1
;  params[fixed_params ne 0] are the values the parameters are fixed at
; image= image to fit
;
;RETURNS:
; set of initial parameters that make sense for the image
;____________________________________






FUNCTION init_cond, params, fixed_params, image, fracs=fracs, rescale=rescale, $
                    free_coords=free_coords

  plist = indgen(n_elements(params)/8)*8
  nprof = n_elements(plist)
  i0_ind = 0
  r_ind = 1
  n_ind = 2
  q_ind = 3
  x0_ind = 5
  phi_ind = 7
  inits = params


  if ((keyword_set(free_coords) eq 0) and $
      (where(fixed_params[plist+x0_ind] eq 1, /null) ne !NULL)) then begin
     setx0 = (where(fixed_params[plist+x0_ind] eq 1))[0]
     inits[plist+x0_ind] = inits[plist[setx0]+x0_ind]
     inits[plist + x0_ind+1] = inits[plist[setx0]+x0_ind+1]
     fixed_params[plist+x0_ind] = 1
     fixed_params[plist+x0_ind+1] = 1
  endif else begin
     if ((keyword_set(free_coords) eq 0) and $
        (where(inits[plist+x0_ind] ne 0, /null) ne !NULL)) then begin
        setx0 = (where(inits[plist+x0_ind] ne 0))[0]
        inits[plist + x0_ind] = inits[plist[setx0] + x0_ind]
        inits[plist + x0_ind+1] = inits[plist[setx0] + x0_ind+1]
     endif else begin
        xs = getXY_start(image)
        resetx0 = where((inits[plist+x0_ind] eq 0) and $
                        (fixed_params[plist+x0_ind] eq 0), /null)
        inits[plist[resetx0]+x0_ind] = xs[0]
        inits[plist[resetx0]+x0_ind+1] = xs[1]
     endelse
  endelse

  setscaling, image, inits, fixed_params, rescale=rescale, fracs=fracs

  return, inits
  

END
