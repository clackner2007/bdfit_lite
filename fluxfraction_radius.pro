;____________________________________
;
; Claire Lackner
;  Jul 14 2010
;
;  calculates the semimajor axis of the ellipse that contains 80% of
;  the flux from the galaxy. The shape of the ellipse is determined by
;  the dominate galaxy component
;
;____________________________________

;returns the difference between a fraction of the flux and the flux
;contained within a given radius ellipse, given the ellipse
FUNCTION flux_diff, rad

COMMON flux_share, enc_flux, im, ellipse

flux=dblarr(n_elements(rad))

for i=0L,n_elements(rad)-1 do begin
    flux[i] = total(im*(ellipse le rad[i]))
end
;print, flux
;print, enc_flux
return, flux - enc_flux

END

;the initial guess is 1.78*r_eff of the dominant component
;the fluxlevel is the flux i want enclosed in the ellipse (80% total)
FUNCTION fluxfraction_radius, image, fluxlevel=fluxlevel, q=q, phi=phi, x0=x0, y0=y0, guess=guess

COMMON flux_share, fl, im, ellp

fl=fluxlevel
im=image

nx = (size(image))[1]
ny = (size(image))[2]

;make ellipse with distance marked
dist_ellipse, ellp, [nx,ny], x0, y0, 1.0/q, phi*180.0/!pi-90.0


;secant method
;result_rad = secant('flux_diff', guess[0], guess[1], precision=0.1)

;bisection
while( flux_diff(guess[0]) gt 0.0 ) do guess[0] *= 0.7
while( flux_diff(guess[1]) lt 0.0 ) do guess[1] *= 1.5
result_rad = bisection('flux_diff', guess[0], guess[1], precision=0.01)

return, result_rad

END

