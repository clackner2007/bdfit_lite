;____________________________________
;
; Claire Lackner
;  Jan 14 2010
;
;  get_imivarpsf.pro
;
; returns a data object with the COSMOS image, ivar, and psf
; uses and exposure time of 2028.0 seconds to get image variance right
;  THERE MAY STILL BE ODD FACTORS OF ~4 hanging around from drizzling
;____________________________________


FUNCTION get_imivarpsf, id, name, rootpath


if file_test(rootpath+'images/'+name+'.fits') eq 0 then begin
    print, "couldn't find file "+rootpath+"images/"+name+'.fits'
    return, {image:0.0, ivar:0.0, psf:0.0, mask:0.0}
 endif

if file_test(rootpath+'ivar/'+name+'.wht.fits') eq 0 then begin
    print, "couldn't find file "+rootpath+"ivar/"+name+'.wht.fits'
    return, {image:0.0, ivar:0.0, psf:0.0, mask:0.0}
 endif

if file_test(rootpath+'masks_all/'+$
             string(id, format='(i0)')+'_mask.fits') eq 0 then begin
    print, "couldn't find file "+$
           rootpath+string(id, format='(i0)')+'_mask.fits for '+name
    return, {image:0.0, ivar:0.0, psf:0.0, mask:0.0}
endif
im = mrdfits(rootpath+'images/'+name+'.fits', 0)
ivar0 = mrdfits(rootpath+'ivar/'+name+'.wht.fits',0)
ivarMask = mrdfits(rootpath+'masks_all/'+string(id, format='(i0)')+'_mask.fits', 0)
mymask = mrdfits(rootpath+"masks_all/"+string(id, format='(i0)')+'_seg.fits', 0)

;temporarily change psf to random psf from list:
;idpath = string(14292, 149.543219D, 2.152322D, $
;                format='(i0,".0_",d10.6,"_",d8.6)')
;print, rootpath+'psf/'+idpath+'.psf.fits.gz'
;locpsf = mrdfits(rootpath+'psf/'+idpath+'.psf_est_red_halo.fits',0)
locpsf = mrdfits(rootpath+'psf/PSF_small.fits',0)
locpsf /= total(locpsf)
mask = mymask
mask[where(mymask-ivarMask ne mymask)] = 0
mask[where(mask ne 0)] = 1

;ivarMask[where(ivarMask gt 0)] = 1
;ivarMask=smooth(ivarMask,50)
;ivar0[where(ivarMask lt 1.e-15)] = 0


;do the image cropping here

;add the variance from the image to variance maps, not done in cosmos
;but everything is in electrons/seconds, and the exposure time is 2028
;sec
; for the ACS images
ivar = ivar0
flip = where((ivar0 ne 0) and (mymask ne 0))
ivar[flip] = 1./(1./ivar0[flip] + im[flip]/2028.0)
ivar[where(ivarMask eq 1)] = 0.0
im[where(ivarMask eq 1)] = 0.0

;mask negative flux, these are noise dominated pixels
nulls=where(ivar lt 0.0,n_nulls)
print, n_nulls
if n_nulls gt 0 then $
  ivar[nulls] *= 0.0

data = {image:im, ivar:ivar, psf:locpsf, mask:mask}

return, data


END



; FUNCTION get_imivarpsf2, id, ra, dec, rootpath

; idpath = string(id, ra, dec, format='(i0,".0_",d10.6,"_",d8.6)')

; ;if file_test(rootpath+'images/'+idpath+'_processed.fits.gz') eq 0 then begin
; ;    print, "couldn't find file "+rootpath+"images/"+idpath
; ;    return, {image:0.0, ivar:0.0, psf:0.0}
; ;endif

; im = mrdfits('RM_COSMOS/0006_skyTest/0006_ring2Image.fits',0)
; ivar0 = mrdfits('RM_COSMOS/0006_skyTest/0006_ring2Ivar.fits',0)
; ;im = mrdfits(rootpath+'images/'+idpath+'_processed.fits.gz',0)
; ;ivarMask = mrdfits(rootpath+'ivar/'+idpath+'.wht.mask.fits',0)
; ;ivar0 = mrdfits(rootpath+'ivar2/'+idpath+'.wht.fits',0)


; locpsf = mrdfits(rootpath+'psf/'+idpath+'.psf.fits.gz',0)


; ;ivarMask[where(ivarMask gt 0)] = 1
; ;ivarMask=smooth(ivarMask,50)
; ;ivar0[where(ivarMask lt 1.e-15)] = 0


; ;do the image cropping here

; ;add the variance from the image to variance maps, not done in cosmos
; ;but everything is in electrons/seconds, and the exposure time is 2028
; ;sec
; ; for the ACS images
; ivar = ivar0
; flip = where(ivar0 ne 0)
; ivar[flip] = 1./(1./ivar0[flip] + im[flip]/2028.0)

; ;mask negative flux, these are noise dominated pixels
; nulls=where(ivar lt 0.0,n_nulls)
; if n_nulls gt 0 then $
;   ivar[nulls] *= 0.0

; data = {image:im, ivar:ivar, psf:locpsf}

; return, data


; END
