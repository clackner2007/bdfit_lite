;+-------------------------------
;NAME:
; model_image
;PURPOSE:
; saves an image of the model fit with a given size (usually same as
; the real image
;
;CALL:
; model_image, sizex, sizey, params, psfImage, cutoff=cutoff
;
;RETURN:
; model_image of sizex x sizey
;
;-------------------------------

FUNCTION model_image, sizex, sizey, params, psfImage, cutoff=cutoff
  
  if not keyword_set(cutoff) then cutoff=0

  makegrid, sizex, sizey, x, y
  model = pixelfluxpsf(x, y, params, _EXTRA={cutoff:cutoff, psfImage:psfImage})

  return, model

END

