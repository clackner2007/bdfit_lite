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


FUNCTION cropimage, image, ivar, psf, mask, x0=x0, y0=y0

nx = (size( image, /dimensions))[0]
ny = (size( image, /dimensions))[1]

makegrid, nx, ny, x, y

mask_y = total(total(mask, 1), /cumulative)/$
  (1.0*total(mask))
mask_x = total(total(mask, 2),/cumulative)/$
  (1.0*total(mask))

x_min = (where(mask_x gt 0.02))[0]
x_max = (where(mask_x gt 0.98))[0]
y_min = (where(mask_y gt 0.02))[0]
y_max = (where(mask_y gt 0.98))[0]


x_mid = (where(mask_x gt 0.5))[0]
y_mid = (where(mask_y gt 0.5))[0]

;print, x_min, y_min, x_max, y_max, x_mid, y_mid

grow=2.0
left = max([0, x_min - fix((x_mid-x_min)*grow)])
low = max([0, y_min - fix((y_mid-y_min)*grow)])
right = min([nx-1, x_max + fix((x_max-x_mid)*grow)])
high = min([ny-1, y_max + fix((y_max-y_mid)*grow)])


;crop based on weights maps
;keep = where(ivar ne 0)
;left = min(x[keep])
;right = max(x[keep])
;low= min(y[keep])
;high = max(y[keep])

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
