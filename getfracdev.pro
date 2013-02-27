;calculate the frac dev from the deVaucouleurs and exponential fits
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;scalefit:
; ARGUMENTS: image, psf, ivar = image, psf, ivar in band to be fit
;            expParams = parameters for exp fit
;            dvcParams = parameters for dvc fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION getfracdev, image, psf, ivar, expParams, dvcParams

imsize = size(image, /dimensions)
makegrid, imsize[0], imsize[1], x,y


disk = pixelfluxpsf(x,y, expParams,_EXTRA={psfImage:psf, cutoff:0})
bulge =pixelfluxpsf(x,y, dvcParams,_EXTRA={psfImage:psf, cutoff:0})

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;linear method
dd=total(ivar*disk*disk)
bb=total(ivar*bulge*bulge)
bd=total(ivar*bulge*disk)

d_norm = (total(ivar*disk*image)*bb - total(ivar*bulge*image)*bd) / $
  (dd*bb-bd^2)
b_norm = (total(ivar*bulge*image)*dd - total(ivar*disk*image)*bd) / $
  (dd*bb-bd^2)

factors = [d_norm, b_norm]


    
;;;;;;QR method;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;flat_isigma = reform(sqrt(ivar),n_elements(ivar))
;    
;design_mat = [[reform(disk,n_elements(disk))*flat_isigma],$
;              [reform(bulge,n_elements(bulge))*flat_isigma]]
;qrfac, design_mat, r, ipvt, /pivot
;qrparams = qrsolv(design_mat, r, $
;                  reform(image, n_elements(image))*flat_isigma, $
;                  pivots=ipvt)
;
;factors = qrparams

; if (factors[1] ne factors[1]) or (factors[0] ne factors[0]) then begin
; ;;;;;;;;;;;   SVD method ;;;;;;;;;;;;;;;;;;;;;;
    
;     flat_isigma = reform(sqrt(ivar),n_elements(ivar))
    
;     design_mat = [[reform(disk,n_elements(disk))*flat_isigma],$
;                   [reform(bulge,n_elements(bulge))*flat_isigma]]

;     svdc, design_mat, w, u, v, /column, /double
;     b = reform(image,n_elements(image))*flat_isigma
;     inv_w = dblarr(2)
;     nonzeros=where(abs(w) < max(w)*1.0e-18*n_elements(b),n)
;     if( n gt 0 ) then $
;       inv_w[nonzeros] = 1.0/w[nonzeros]
;     factors = ((u##b)*inv_w)#v

;     inv_w_mat = [[inv_w^2],[inv_w^2]]
;     covar=(inv_w_mat*v)##v
;endif

;print, factors
if factors[1] lt 0 then return, 0.0
if factors[1] gt 1.0 then return, 1.0
return, factors[1]

;return, factors[1]/(factors[0] + factors[1])

END
