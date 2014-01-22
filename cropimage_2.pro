;____________________________________
;
; Claire Lackner
;  Sept 17, 2009
;
;  cropimage.pro
;
; crops image such that only relvant object is in filed, does so using
; a flux threshold
;
;____________________________________


FUNCTION cropimage_2, image, ivar, psf, mask, x0=x0, y0=y0

nx = (size( image, /dimensions))[0]
ny = (size( image, /dimensions))[1]

makegrid, nx, ny, x, y

good = where(mask gt 0)
if (size(good, /dimensions) eq 0) then begin
   return, {image:0.0, ivar:0.0, psf:0.0, mask:0.0}
endif
totgood = n_elements(good)
left = min(good mod nx)
right = max(good mod nx)
low = (good / nx)[0]
high = (good / nx)[totgood-1]


totsky = n_elements(where((mask[left:right,low:high] eq 0) and $
                          (ivar[left:right,low:high] ne 0.0)))
ratiosky = 2
adjust = 0

while ((totsky lt ratiosky*totgood) and $
       ((left gt 0) or (right lt nx-1) or (low gt 0) or $
        (high lt ny-1))) do begin
   print, totsky, totgood*ratiosky, left, right, low, high
   xmid = (right - left)/2 + left
   ymid = (high - low)/2 + low
   q = (high - low)*1.0 / (right-left)
   totarr = (ratiosky +adjust*0.1 + 1)*totgood
  
   left = max([min([left-1, xmid - fix(sqrt(totarr/q)*0.5)]), 0])
   right = min([max([right+1, xmid + fix(sqrt(totarr/q)*0.5)]), nx-1])
   low = max([0, min([low-1, ymid - fix(sqrt(totarr*q)*0.5)])])
   high = min([ny-1, max([high+1,ymid + fix(sqrt(totarr*q)*0.5)])])
   totsky = n_elements(where((mask[left:right,low:high] eq 0) and $
                             (ivar[left:right,low:high] ne 0.0)))
   adjust += 1
endwhile


x0 = left
y0 = low

;set a minimun image size
if ((left ge right) or (low ge high)) then $
   return, {image:0.0, ivar:0.0, psf:0.0, mask:0.0}

cutout = dblarr(right-left+1, high-low+1)
cutout = image[left:right, low:high]

vcutout = dblarr(right-left+1,high-low+1)
vcutout = ivar[left:right, low:high]
;psfcutout = cutout*0.00

maskcutout = dblarr(right-left+1, high-low+1)
maskcutout = mask[left:right, low:high]

return, {image:cutout, ivar:vcutout, psf:psf, mask:maskcutout}

END
