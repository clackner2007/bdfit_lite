;makes a one d image convolved with the psf

FUNCTION pixelfluxpsf1d, x, params, _EXTRA=_EXTRA

makegrid, n_elements(x)*2, n_elements(x)*2, xx, yy

image = pixelfluxpsf(xx,yy,params,_EXTRA=_EXTRA)

oned_image, flux, rad, ivar, image, _EXTRA.pa, _EXTRA.q, $
  xcen=params[12],ycen=params[13]

return, flux


END
