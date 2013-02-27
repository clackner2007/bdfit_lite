;____________________________________
;
; Claire Lackner
;  Sept 7, 2009
;
;  makepsf.pro
;
;  makes an image of the psf, given the fwhm and using the image
;  dimsions given, this also normalizes the psf such that the total
;  flux is 1.0
;   
;____________________________________

FUNCTION makepsf, nx, ny, fwhm

c = dindgen(nx)
r = dindgen( ny )
;set up x and y at integer values
x = c#(r*0+1)
y = (c*0+1)#r
;center of gaussian and image
xcen = nx/2.0D
ycen = ny/2.0D
;sigma, related to fwhm
sigma = fwhm / (2.*sqrt(2.*alog(2.0)))

psf = gauss2( x, y, [xcen, ycen, sigma, 1.0] )
psf /= total(psf) ;normalize total flux to 1

;thresh = where( psf lt 1.0e-15, count )

;if count gt 0 then psf[thresh] = 0.0D

;psf /= total(psf)

return, psf
END
