;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get the chi^2 value for a given model for a given image
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION getChiSquare, image, ivar, model_params, _EXTRA=_EXTRA

makegrid, size(image, /dimension)[0], size(image, /dimension)[1], x, y

modelimage = pixelfluxpsf(x,y,model_params,_EXTRA=_EXTRA)

return, total((image-modelimage)**2*ivar)

END


