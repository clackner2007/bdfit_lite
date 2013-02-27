;Claire Lackner
; Jul 16 2010
;
; get_r80.pro
;
; function to return the radius containing 80% of the flux of the
; galaxy
;
; the method uses the non-cutoff value of the flux such that the 1/2
; light radius is actually the half light radius, using the cutoff
; value lowers the answer by about 1% (<0.5 pixels for many galaxies)
;

;function for best guess, this is given by linear formula for r_80 of
;single sersic component
FUNCTION get_guess, n, reff
;return, (n*0.46804314 + 1.31185517)*reff
return, (0.002819*n^3 - 0.02982*n^2 + 0.5567*n + 1.248)*reff
END

FUNCTION get_r80, fit_params

if( n_elements(fit_params) le 8 ) then $
  return, get_guess(fit_params[2], fit_params[1])

q=fit_params[3]
phi=fit_params[7]
guess=get_guess(fit_params[2],fit_params[1])
dimen = round(1.7*guess)*2

if(n_elements(fit_params) gt 8) then begin
    if(bulgetotot(fit_params) gt 0.5 ) then begin
        q=fit_params[9]
        phi=fit_params[11]
        guess=get_guess(fit_params[8],fit_params[7])
    endif else phi = fit_params[5]
    dimen = round(1.7*guess)*2
    fit_params[12] = dimen/2
    fit_params[13] = dimen/2
endif else begin
    fit_params[5] = dimen/2
    fit_params[6] = dimen/2
endelse

makegrid, dimen,dimen,x,y

model = pixelfluxpsf(x,y,fit_params,_EXTRA={cutoff:0})
flux80 = 0.8*(totalsersicflux([fit_params[0:4],0.0,0.0,0.0],cutoff=0) + $
              totalsersicflux([fit_params[6:10],0.0,0.0,0.0],cutoff=0))
r80=fluxfraction_radius(model, fluxlevel=flux80, $
                            q=q,phi=phi,x0=dimen/2,y0=dimen/2,$
                            guess=[min([fit_params[1],fit_params[7]]),guess])

;print, guess, r80, abs(guess-r80)/r80
return, r80
END
