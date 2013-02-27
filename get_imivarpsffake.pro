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


FUNCTION get_imivarpsffake, id, ra, dec, rootpath, isnr

idpath = id ;, ra, dec, format='(i0,".0_",d10.6,"_",d8.6)')

if file_test(rootpath+'images/'+idpath) eq 0 then begin
    print, "couldn't find file "+rootpath+"images/"+idpath
    return, {image:0.0, ivar:0.0, psf:0.0}
endif
im = mrdfits(rootpath+'images/'+idpath,0)
locpsf = mrdfits(rootpath+'SNR_inf.psf.fits',0)

ivar=dblarr((size(im))[1],(size(im))[2]) + 1.0
if isnr ne 0 then $
  ivar *= isnr
mask=ivar
size_im = size(im, /dimensions)
;take the inner 1/5, crop image will double size
;factor=11
;oneT = size_im[0]/factor
;ff = factor/2
;print, factor, ff, oneT
;mask[0:ff*oneT,*] = 0.0
;mask[((factor-ff)*oneT):size_im[0]-1,*]=0.0
;mask[ff*oneT:(factor-ff)*oneT,0:ff*oneT]=0.0
;mask[ff*oneT:(factor-ff)*oneT,((factor-ff)*oneT):size_im[1]-1]=0.0

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
