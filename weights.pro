;____________________________________
;
; Claire Lackner
;  Sept 17, 2009
;
;  weights.pro
;
;  IDL
;    calculates the set of weights to be used for each pixel in
;    fitting a sersic profile, this overrides the error image given by
;    sdss, although that image can still be used
;____________________________________


FUNCTION weights, ivar, image, params, radius

nx = (size( image, /dimensions))[0]
ny = (size( image, /dimensions))[1]

w = dblarr( nx, ny )

;ivar weights
;w += ivar

;weights based on assumed fit
;  row = dindgen( nx )
;  col = dindgen( ny )
;  x = row#(col*0+1)
;  y = (row*0+1)#col

;  w += ivar*sersic( x, y, params )/total(sersic(x,y,params))

;weights such that only a smaller, central region is fit
;if the radius is set to zero, the entire image is fit
;weights zero outside central region
if( radius gt 1e-5 ) then begin
    w = make_circ( w, radius+1., cent=[params[5], params[6]] )*ivar
endif else begin
    w += ivar
endelse

return, w

END
