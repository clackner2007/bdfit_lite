;+
;NAME:
;   sdssize_image
;PURPOSE:
; takes a cosmos image and returns an image
; more simlar to sdss image (don't worry about details of tiny tim
; psf)
;
;CALLING SEQUENCE:
;   sdssize_image( image, psf, ivar, newsigma, rebinsize, added_noise,
;                  orig_psf_sigma= )
;
;INPUTS:
; image_orig ==original COSMOS image
; psf_orig == original COSMOS psf
; ivar_orig == original COSMOS ivar map (including galaxy?)
; mask_orig == masked ivar of image
; new_psf_sigma == width of new psf (in NEW image pixels)
; add_noise_factor == factor by which to decrease the signal-to=noise
; (always greater than 1) 
; rebinsize == number of pixels to rebin (always greater than 1)
; orig_psf_sigma == COSMOS psf width in COSMOS pixels (default=1.75)
;
;OUTPUTS:
; data structure with new image, ivar and psf
;
;HISTORY:
; Claire Lackner
; March 2012
;
;-


FUNCTION sdssize_image, image_orig, psf_orig, ivar_orig, $
                        mask_orig, $
                        new_psf_sigma, rebinsize, add_noise_factor, $
                        orig_psf_sigma=orig_psf_sigma


orig_size = size(image_orig, /dimensions)

;do things to the var map instead of the ivar map
var_orig = dblarr(orig_size[0], orig_size[1])
var_orig[where(ivar_orig ne 0)] = $
  1.0/ivar_orig[where(ivar_orig ne 0)]



;first convolve with larger psf
if not keyword_set(orig_psf_sigma) $
  then orig_psf_sigma=1.75

;in the old units, the size of the psf to be added is:
npsfsig=sqrt((new_psf_sigma*rebinsize)^2 - orig_psf_sigma^2)

add_psf = psf_gaussian( npixel=101, fwhm=npsfsig*2.35, $
                        ndimen=2, /double, /normalize )

new_image = convolve(image_orig, add_psf, $
                     ft_psf=ft_add_psf, $
                     /correlate)
new_psf = convolve(psf_orig, add_psf, $
                   ft_psf=ft_add_psf, $
                   /correlate)

;this assumes that the variances are uncorrelated which is wrong
; but will have to do for now
new_var = var_orig /(4.0*!pi*npsfsig^2)

;;;now, rebin the pixels
new_image = congrid( new_image, orig_size[0]/(rebinsize*1.0),$
                     orig_size[1]/(rebinsize*1.0), $
                     /center, /interp ) * rebinsize^2

psf_size = size(new_psf, /dimensions)
new_psf = congrid( new_psf, psf_size[0]/(rebinsize*1.0),$
                     psf_size[1]/(rebinsize*1.0), $
                     /center, /interp ) * rebinsize^2

;because the flux is conserved...the variance stays the same
;again this requires that the variance is not correlated
new_var = congrid( new_var, orig_size[0]/(rebinsize*1.0),$
                   orig_size[1]/(rebinsize*1.0), $
                   /center, /interp ) * rebinsize^2

;bin up the mask, but don't do anything else to it
new_mask = congrid( mask_orig, orig_size[0]/(rebinsize*1.0),$
                   orig_size[1]/(rebinsize*1.0), $
                   /center, /interp ) * rebinsize^2

;the original variance in the ivar map is too large compared
; to the actual image


;finally...deal with the noise
var_add = new_var * (add_noise_factor^2 * $
                     4.0 * !pi * npsfsig^2 - 1 )
pp=where(var_add lt 0, badvar)
if( badvar gt 0) then $
  print, "some of the variances are less than zero"

new_var = new_var + var_add
new_ivar = new_var
new_ivar[where(new_ivar ne 0)] = $
  1.0/new_ivar[where(new_ivar ne 0)]

new_size = size(new_image, /dimensions)
noise_add = randomn(seed, new_size[0], new_size[1]) * $
  sqrt(var_add)

new_image += noise_add


;return the whole kit a kaboodle
return, {image:new_image, ivar:new_ivar, $
         psf:new_psf, mask:new_mask}

END
