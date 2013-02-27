;____________________________________
;
; Claire Lackner
;  Sept 7, 2009
;
;  pixelflux.pro
;
;  IDL
; takes the model+psf convolution and returns the value in a given
; pixel, at a position (px, py) 
;      
; if the gradient in a pixel is above a certain threshold, the pixel,
; is subdivided into four mini-pixels, which continues until the
; gradient of the flux is below a threshold and the mini-pixel fluxes
; are summed together 
;
;  the pixel is indexed by integers in the center, and has a size of
;  1.0x1.0 
;   
;____________________________________


;function for the flux, assuming there are multiple sersic profiles
;assuming they are all centered at the same point, then, the position
;parameters are at the end
FUNCTION profileFlux, px, py, params, cutoff=cutoff, rescale=rescale

if size(px, /n_dimensions) gt 1 then begin
    flux = dblarr( (size( px, /dimensions))[0], $
                   (size( py, /dimensions))[1] )
endif else begin
    flux = 0.0D
endelse

if( n_elements(params) eq 8 ) then flux = sersic( px, py, params, cutoff=cutoff ) $
;else begin
;    numprof = (n_elements(params)-2)/6
;    last= n_elements(params)
;    for i=0,numprof-1 do begin
;        p = dblarr(8)
;        p[0:4] = params[i*6:i*6+4]
;        p[5:6] = params[last-2:last-1]
;        p[7] = params[i*6+5]
;        flux += sersic( px, py, p, cutoff=cutoff )
;    endfor
else begin
    numprof=(n_elements(params))/8
    last = n_elements(params)
    for i=0,numprof-1 do begin
        p=dblarr(8)
        p = params[i*8:i*8+7]
        if keyword_set(rescale) and (i gt 0) then  p[1] *= params[1]
        flux += sersic(px, py, p, cutoff=cutoff)
    endfor
endelse

return, flux

END


;function that computes the flux of the sersic model in a given pixel
;assuming that the psf is a delta function
; this is the function that is passed through to mpfitfun, for a delta
; function psf
FUNCTION pixelflux, px, py, params, cutoff=cutoff, imax=imax, rescale=rescale


if not keyword_set(imax) then imax=8

THRESH = 1.0e-2

nx = (size( px, /dimensions))[0]
ny = (size( py, /dimensions))[1]
flux0 =dblarr(nx,ny)
flux0 = profileFlux( px, py, params, cutoff=cutoff, rescale=rescale )
x0 = px[0]
y0 = py[0]

i=1

bad = [1]
while i lt imax and bad[0] ne -1 do begin
;    print, n_elements(bad)
    row = (0.5+dindgen(nx*2^i))/(2^i) + x0 - 0.5
    col = (0.5+dindgen(ny*2^i))/(2^i) + y0 - 0.5
    
    x = row#(col*0+1.)
    y = (row*0+1)#col

    
    flux_large = profileFlux( x, y, params, cutoff=cutoff, rescale=rescale )

    fluxbin = rebin( flux_large, nx, ny )

    ;put the binned flux into the old flux 
    fluxnew = flux0
    fluxnew[x0-px[0],y0-py[0]] = fluxbin

    
    ;take the difference and selecte the pixels that didn't converge
    ;the 2.5e-3 comes from ignoring the edges where the flux starts to
    ;get seriously truncated, which is to speed up the calculation
    trunc_diff=1.e-8
    if keyword_set(cutoff) then trunc_diff=2.8e-3
    if( max(flux0) lt 1.0e-12 ) then return, flux0
    notzeros = where( abs(flux0/max(flux0)) gt trunc_diff, count )
    difflux = abs( flux0 - fluxnew )
    if( count gt 0 ) then difflux[notzeros] = difflux[notzeros]/flux0[notzeros]
    bad = where((difflux gt THRESH) and (abs(flux0/max(flux0)) gt trunc_diff))
    
    ;get out of while loop if there are no elements under the 
    ;threshold
    if bad[0] eq -1 then begin
;     print, 'getting out'
        continue
    endif
   ;set the new parameters for the area that need more subdividing
    x0 = min( px[bad] )
    y0 = min( py[bad] )
    nx = max( px[bad] ) - x0 + 1
    ny = max( py[bad] ) - y0 + 1
    
    ;copy fluxnew into fluxold
    flux0 = fluxnew
    
    i++
endwhile

;if( i ge 8 and bad[0] ne -1) then $
;  print, 'flux did not converge',params[0:5]
;if( i ge 8 and bad[0] ne -1 and n_elements(params) gt 13) then $
;  print, " ",params[6:11]

return, flux0

END

FUNCTION pixelfluxpsf, x, y, params, imax=imax, _EXTRA=_EXTRA

nparams = n_elements(params)
;check there is a sky value,
if nparams mod 8 eq 0 then begin
    skyVal = 0.0
    sparams = params
endif else begin
    skyVal = params[nparams-1]
    sparams = [params[0:nparams-2]]
endelse

rescale=0
if tag_exist(_EXTRA, 'rescale', /top_level) then rescale=_EXTRA.rescale

if( not tag_exist(_EXTRA, 'psfImage', /top_level) ) then begin
  if( tag_exist(_EXTRA, 'cutoff', /top_level ) ) then $
  flux=pixelflux(x,y,sparams,cutoff=_EXTRA.cutoff, rescale=rescale) else $
  flux = pixelflux(x,y,sparams,cutoff=0, rescale=rescale)
  return, flux + skyVal
endif

psfcrop = _EXTRA.psfImage

                                        ;size of oringal psf
 xpsf = (size( psfcrop, /dimensions ))[0]
 ypsf = (size( psfcrop, /dimensions ))[1]
 nx_im = (size(x, /dimensions))[0]
 ny_im = (size(y, /dimensions))[1]
 nx = max([xpsf, nx_im])
 ny = max([ypsf, ny_im])

 psf = dblarr( nx, ny )

 psf[nx-xpsf/2, ny-ypsf/2] = psfcrop[0:(xpsf/2-1), 0:(ypsf/2-1)]
 psf[0, ny-ypsf/2] = psfcrop[xpsf/2:xpsf-1, 0:(ypsf/2-1)]
 psf[nx-xpsf/2, 0] = psfcrop[0:(xpsf/2-1), ypsf/2:ypsf-1]
 psf[0, 0] = psfcrop[xpsf/2:xpsf-1, ypsf/2:ypsf-1]


 flux = dblarr( nx, ny )
 if( tag_exist(_EXTRA, 'cutoff', /top_level ) ) then $
   flux[(nx-nx_im)/2, (ny-ny_im)/2]=pixelflux(x,y,sparams,$
                                              cutoff=_EXTRA.cutoff, $
                                             rescale=rescale) else $
   flux[(nx-nx_im)/2, (ny-ny_im)/2] = pixelflux(x,y,sparams,cutoff=0, rescale=rescale)

 ;flux=reverse(flux,1)
 ;flux=reverse(flux,2)
;  ;test whether simply smoothing the function works, instead of convolving
;  ; it does, suggesting the derivative may not be the biggest problem
;  ;flux = smooth( flux, [4,4] )

;  ;fft
 im_fft = fft( flux, /double)
 psf_fft = fft(psf, /double)
 conv_fft = im_fft * psf_fft

 flux = (double(fft(conv_fft,/double, /inverse)))*nx*ny
 ;flux=reverse(flux,1)
 ;flux=reverse(flux,2)
;flux = convolve(flux,psfcrop);,/correlate)

return, flux[(nx-nx_im)/2:(nx-nx_im)/2+nx_im-1, $
             (ny-ny_im)/2:(ny-ny_im)/2 + ny_im-1] + skyVal



END


