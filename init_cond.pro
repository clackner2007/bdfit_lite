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


FUNCTION init_cond, params, fixed_params, image, fracs=fracs, rescale=rescale

  plist = indgen(n_elements(params)/8)*8
  nprof = n_elements(plist)
  i0_ind = 0
  r_ind = 1
  n_ind = 2
  q_ind = 3
  x0_ind = 5
  phi_ind = 7
  inits = params
  rescale keyword_set(rescale)

  resetX = where(fixed_params[plist+x0_ind] eq 0 )
  resetY = where(fixed_params[plist+x0_ind+1] eq 0 )

  if (resetX ne 0) or (resetY ne 0) then begin
     xlen = (size(image, /dimensions))[0]
     ylen = (size(image, /dimensions))[1]
     ys = reverse(sort(smooth(total(data.image,1),min([4, ylen-1]))))
     xs=reverse(sort(smooth(total(data.image,2),min([4, xlen-1]))))
     y0=ys[0]
     x0=xs[0]
     inits[plist[resetX]+x0_ind] = x0
     inits[plist[resetY]+y0_ind] = y0
  endif

  tot = max([total(data.image),10])
  totp = dblarr(n_elements(plist))
  for ip=0, n_elements(plist)-1 do begin
     totp[ip] = max(totalsersicflux(inits[plist[ip]:plist[ip]+7]), 1.0e-3)
  endfor

  if n_elements(fracs) eq 0 then fracs = dblarr(nprof)+1.0/nprof
  for ip=0, n_elements(plist)-1 do begin
     if fixed_params[plist[ip]+i0_ind] ne 0 then continue
     inits[plist[ip]] = frac[ip]*tot/totp[ip]
  endfor
  
  
  return, inits
  

END

;if(keyword_set(nodisk)) then parinfo[*].value = [1.0,10.,bulgeSersic, $
;                                                0.7,0.0,x0,y0,0.2, 0.0] $
;else begin
;    parinfo[*].value = [1.0, 30.1, diskSersic, 0.85, $
 ;                       0.0, x0, y0, 0.2, $
 ;                       1.0, 8.1, bulgeSersic, 0.7, 0.0, $
 ;                       x0, y0, 0.2, 0.0]
    

frac=0.5

if( keyword_set(_EXTRA) ) then begin

    if tag_exist(_EXTRA, 'start') then parinfo[*].value=_EXTRA.start $
    else begin
        
        if(tag_exist(_EXTRA,"Reff")) then begin
            parinfo[1].value = _EXTRA.Reff
            if(not keyword_set(nodisk)) then begin
               if rescale eq 0 then $
                  parinfo[9].value = _EXTRA.Reff*0.5 $
               else parinfo[9].value = 0.5
                frac=0.5
                if tag_exist(_EXTRA,"fracdev") then begin
                    frac=(0.10)+_EXTRA.fracdev*0.8
                    if(_EXTRA.fracdev lt 0.5) then $
                      parinfo[9].value *= 0.5 $
                    else $
                      parinfo[1].value *= 1.0 ;0.5 + _EXTRA.fracdev ;-1.+2.*_EXTRA.fracdev
                endif

            endif 
        endif 
        if(tag_exist(_EXTRA,"q")) then begin
            parinfo[3].value = max([0.11,_EXTRA.q])
            if(not keyword_set(nodisk)) then parinfo[11].value = max([0.11,_EXTRA.q])
        endif
        if(tag_exist(_EXTRA,"phi")) then begin
            parinfo[phi_ind].value = _EXTRA.phi
            if(not keyword_set(nodisk)) then parinfo[15].value = _EXTRA.phi
        endif
    endelse
endif 
                                        ;to maximum pixel flux
tot = max([total(data.image),10])
if( nparams lt 16) then $
  parinfo[0].value = tot/totalsersicflux(start_params) $
else begin
   rb=start_params[9]
   if rescale ne 0 then rb *= start_params[1]
    d = max([totalsersicflux([start_params[0:4],0.0,0.0,0.0]),1.0e-2])
    b = max([totalsersicflux([start_params[8], $
                              rb, $
                              start_params[10:12],0.0,0.0,0.0]),1.0e-2])

    parinfo[8].value = frac*tot / b
    parinfo[0].value = tot*(1.0-frac)/d
endelse
  start_params = parinfo[*].value
 
