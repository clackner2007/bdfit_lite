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
;
; Assumes the data is arranged as:
; 
;  image = rootpath/images/{id}.0_{ra}_{dec}.processed.fits.gz
;  inverse variance (weights) = rootpath/ivar/{id}.0_{ra}_{dec}.wht.fits
;  mask = rootpath/mask/{id}.0_{ra}_{dec}.wht.mask.fits
;  psf = rootpath/psf/{id}.0_{ra}_{dec}.psf.fits.gz
; 
;
;
;____________________________________


FUNCTION get_imivarpsf, id, ra, dec, rootpath

idpath = string(id, ra, dec, format='(i0,".0_",d10.6,"_",d8.6)')
extra=''

if file_test(rootpath+'images/'+idpath+extra+'_processed.fits.gz') eq 0 then begin
    print, "couldn't find file "+rootpath+"images/"+idpath+extra+'_processed.fits.gz'
    return, {image:0.0, ivar:0.0, psf:0.0}
endif
im = mrdfits(rootpath+'images/'+idpath+extra+'_processed.fits.gz',0)
ivarMask = mrdfits(rootpath+'mask/'+idpath+extra+'.wht.mask.fits', 0);.gz',0)
ivar0 = mrdfits(rootpath+'ivar/'+idpath+extra+'.wht.fits',0);.gz',0)

locpsf = mrdfits(rootpath+'psf/'+idpath+extra+'.psf.fits.gz',0)
mask = ivarMask
mask[where(mask ne 0)] /= mask[where(mask ne 0)]


ivar = ivar0
flip = where(ivar0 ne 0)
ivar[flip] = 1./(1./ivar0[flip] + im[flip]/2028.0)

;mask negative flux, these are noise dominated pixels
nulls=where(ivar lt 0.0,n_nulls)
if n_nulls gt 0 then $
  ivar[nulls] *= 0.0

data = {image:im, ivar:ivar, psf:locpsf, mask:mask}

return, data


END

