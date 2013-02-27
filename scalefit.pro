;linear scale the multiplicative coefficients to fit an image
;the fidicual fit is given by params
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;scalefit:
; ARGUMENTS: image, psf, ivar = image, psf, ivar in band to be fit
;            fitparams = parameters of fit to be rescaled
;                    Example: to fit g using r band fit
;                    fitcolors, image_g, psf_g, ivar_g, params_r
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION scalefit, image, psf, ivar, params, covar=covar

imsize = size(image, /dimensions)
makegrid, imsize[0], imsize[1], x,y

if(n_elements(params) eq 8 ) then begin
    func = pixelfluxpsf(x,y,params,_EXTRA={psfImage:psf,cutoff:1})
    new_norm = total(image*ivar*func) / total(ivar*func*func)
    new_norm *=  params[0]
;    print,  'linear chi_sq ', $
 ;     total(ivar*(image-pixelfluxpsf(x,y,[new_norm,params[1:7]],$
 ;                                    _EXTRA={psfImage:psf,cutoff:1}))^2)

    covar = 1./total(ivar*func*func)
    return, [new_norm,params[1:7]]
    
endif else begin
    bulgeparams = [params[6:10], params[12:13], params[11]]
    diskparams = [params[0:4], params[12:13], params[5]]

    disk = pixelfluxpsf(x,y,diskparams,_EXTRA={psfImage:psf, cutoff:1})
    bulge =pixelfluxpsf(x,y,bulgeparams,_EXTRA={psfImage:psf, cutoff:1})

;;;;;;;;same method, using qr decomposition
    flat_isigma = reform(sqrt(ivar),n_elements(ivar))
    
    design_mat = [[reform(disk,n_elements(disk))*flat_isigma],$
                  [reform(bulge,n_elements(bulge))*flat_isigma]]
    
;     qrfac, design_mat, r, ipvt, /pivot
;     qrparams = qrsolv( design_mat, r, $
;                        reform(image,n_elements(image))*flat_isigma, $
;                        pivots=ipvt )
;     d_norm = qrparams[0]*diskparams[0]
;     b_norm = qrparams[1]*bulgeparams[0]
;    print, 'lin', qrparams*[params[0],params[6]]
;;;;;;;;;;;;;;;;;;;;;;;;;SVD method needs to be used????
;       dd = total(ivar*disk*disk)
;       bb = total(ivar*bulge*bulge)
;       bd = total(ivar*bulge*disk)

;       d_norm = (total(ivar*disk*image)*bb-total(ivar*bulge*image)*bd) /$
;         (dd*bb - bd*bd)
;       b_norm = (total(ivar*bulge*image)*dd-total(ivar*disk*image)*bd) /$
;         (dd*bb-bd*bd)
;;;;;;;;;;;   SVD method ;;;;;;;;;;;;;;;;;;;;;;
     svdc, design_mat, w, u, v, /column, /double
     b = reform(image,n_elements(image))*flat_isigma
     inv_w = dblarr(2)
     nonzeros=where(abs(w) < max(w)*1.0e-18*n_elements(b),n)
     if( n gt 0 ) then $
       inv_w[nonzeros] = 1.0/w[nonzeros]
     factors = ((u##b)*inv_w)#v

     inv_w_mat = [[inv_w^2],[inv_w^2]]
     covar=(inv_w_mat*v)##v
     
     d_norm = factors[0]
     b_norm = factors[1]
     
     d_norm *= diskparams[0]
     b_norm *= bulgeparams[0]

     new_params = [d_norm, params[1:5], b_norm, params[7:13]]
;      print, 'chi-squared ', $
;        total(ivar*(image-pixelfluxpsf(x,y,new_params, $
;                                       _EXTRA={psfImage:psf,cutoff:1}))^2)
     return, new_params

endelse
END
