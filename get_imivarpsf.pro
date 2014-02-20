;____________________________________
;
; Claire Lackner
;  Jan 14 2010
;
;  get_imivarpsf.pro
;
; returns a data object with the image, ivar, psf and a mask (not
;really used)
; this is specfic to the NSA galaxy image format
;
; ARGUMENTS:
; name = IAU name
; rootpath = folder where images live
; atlast id = atlas id from NSA (it's part of filenames)
; parentid = parent object id from NSA (again for filenames)
; band=number [0,1,2,3...] of filter/HDU to use.
;____________________________________


FUNCTION get_imivarpsf, name, rootpath, atlasid, parentid, band=band

;default to r-band image
if n_elements(band) ne 1 then band=2
bands=['u', 'g', 'r', 'i', 'z', 'nd', 'fd']

pstring = strtrim(parentid, 2)


imgname = rootpath+'images/'+name+'-'+pstring + $
          '-atlas-'+strtrim(atlasid, 2)+'.fits.gz' 
if file_test( imgname ) eq 0 then begin
    print, "couldn't find file "+imgname
    return, {image:0.0, ivar:0.0, psf:0.0, mask:0.0}
 endif

ivarname = rootpath+'ivar/'+name+'-parent-'+pstring+'.fits.gz'
if file_test(ivarname) eq 0 then begin
    print, "couldn't find file "+ivarname
    return, {image:0.0, ivar:0.0, psf:0.0, mask:0.0}
 endif

psfname = rootpath+'psf/'+name+'-'+bands[band]+'-bpsf.fits.gz'
if file_test(psfname) eq 0 then begin
    print, "couldn't find file "+psfname
    return, {image:0.0, ivar:0.0, psf:0.0, mask:0.0}
 endif

im = mrdfits(imgname, band)
ivar = mrdfits(ivarname, band*2+1)


;there are a few ways to get the psf, from sdss TODO (add this)
;or from the nsa files
locpsf = mrdfits(psfname, 0)
;normalize psf
locpsf /= total(locpsf)

;make a mask, it's one for the masked out regions
mask = intarr((size(im,/dimensions))[0], $
              (size(im, /dimension))[1])
mask[where(im ne 0)] = 1
ivar[where(mask eq 0)] = 0

data = {image:im, ivar:ivar, psf:locpsf, mask:mask}

return, data

END

