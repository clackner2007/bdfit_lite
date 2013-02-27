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
FUNCTION profileFlux, px, py, params, cutoff=cutoff

if size(px, /n_dimensions) gt 1 then begin
    flux = dblarr( (size( px, /dimensions))[0], $
                   (size( py, /dimensions))[1] )
endif else begin
    flux = 0.0D
endelse

if( n_elements(params) eq 8 ) then flux = sersic( px, py, params, cutoff=cutoff ) $
else begin
    numprof = (n_elements(params)-2)/6
    last= n_elements(params)
    for i=0,numprof-1 do begin
        p = dblarr(8)
        p[0:4] = params[i*6:i*6+4]
        p[5:6] = params[last-2:last-1]
        p[7] = params[i*6+5]
        flux += sersic( px, py, p, cutoff=cutoff )
    endfor
endelse


return, flux

END


;function that computes the flux of the sersic model in a given pixel
;assuming that the psf is a delta function
; this is the function that is passed through to mpfitfun, for a delta
; function psf
FUNCTION pixelflux, px, py, params, cutoff=cutoff


THRESH = 1.0e-4

nx = (size( px, /dimensions))[0]
ny = (size( py, /dimensions))[1]
flux0 =dblarr(nx,ny)
flux0 = profileFlux( px, py, params, cutoff=cutoff )
x0 = px[0]
y0 = py[0]

xstep = px[1]-px[0]
ystep = py[nx]-py[0]

i=1

bad = [1]
while i lt 8 and bad[0] ne -1 do begin
;    print, n_elements(bad)
;    row = (0.5+dindgen(nx*2^i))/(2^i) + x0 - 0.5
;    col = (0.5+dindgen(ny*2^i))/(2^i) + y0 - 0.5
    
;    x = row#(col*0+1.)
;    y = (row*0+1)#col

    makegrid, nx*2^i, ny*2^i, x, y, x_size=xstep/(2^i),$
      y_size=ystep/(2^i), x_start=x0-xstep*0.5*(1.-1./2^i),$
      y_start=y0-ystep*0.5*(1-1./2^i)

    flux_large = profileFlux( x, y, params, cutoff=cutoff )

    fluxbin = rebin( flux_large, nx, ny )

    ;put the binned flux into the old flux 
    fluxnew = flux0
    fluxnew[(x0-px[0])/xstep,(y0-py[0])/ystep] = fluxbin

    
    ;take the difference and selecte the pixels that didn't converge
    ;the 2.5e-3 comes from ignoring the edges where the flux starts to
    ;get seriously truncated, which is to speed up the calculation
    if( max(flux0) lt 1.0e-12 ) then return, flux0
    notzeros = where( abs(flux0/max(flux0)) gt 2.8e-3, count )
    difflux = abs( flux0 - fluxnew )
    if( count gt 0 ) then difflux[notzeros] = difflux[notzeros]/flux0[notzeros]
    bad = where(difflux gt THRESH)
    
    ;get out of while loop if there are no elements under the 
    ;threshold
    if bad[0] eq -1 then begin
;     print, 'getting out'
        continue
    endif
   ;set the new parameters for the area that need more subdividing
    x0 = min( px[bad] )
    y0 = min( py[bad] )
    nx = (max( px[bad] ) - x0)/xstep + 1
    ny = (max( py[bad] ) - y0)/ystep + 1
    
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

FUNCTION pixelfluxpsf, x, y, params, _EXTRA=_EXTRA


if( not tag_exist(_EXTRA, 'psfImage', /top_level) ) then begin
  if( tag_exist(_EXTRA, 'cutoff', /top_level ) ) then $
  flux=pixelflux(x,y,params,cutoff=_EXTRA.cutoff) else $
  flux = pixelflux(x,y,params,cutoff=0)
  return, flux
endif

psfcrop = _EXTRA.psfImage

;factor to oversample before convolving with psf
ff=1
                                        ;size of oringal psf
xpsf = (size( psfcrop, /dimensions ))[0]*ff
ypsf = (size( psfcrop, /dimensions ))[1]*ff
psfcrop = rebin(psfcrop, xpsf, ypsf)/(ff^2)

nx_im = (size(x, /dimensions))[0]*ff
ny_im = (size(y, /dimensions))[1]*ff
nx = max([xpsf, nx_im])
ny = max([ypsf, ny_im])

psf = dblarr( nx, ny )

;padded psf, with the 'negtive' frequencies placed in the opposite edges
psf[nx-xpsf/2, ny-ypsf/2] = psfcrop[0:(xpsf/2-1), 0:(ypsf/2-1)]
psf[0, ny-ypsf/2] = psfcrop[xpsf/2:xpsf-1, 0:(ypsf/2-1)]
psf[nx-xpsf/2, 0] = psfcrop[0:(xpsf/2-1), ypsf/2:ypsf-1]
psf[0, 0] = psfcrop[xpsf/2:xpsf-1, ypsf/2:ypsf-1]

flux = dblarr( nx, ny )
makegrid, nx, ny, xx, yy, x_size = 1./ff, y_size=1./ff, $
  x_start=(-0.5+0.5/ff), y_start=(-0.5+0.5/ff)

if( tag_exist(_EXTRA, 'cutoff', /top_level ) ) then $
  flux[(nx-nx_im)/2, (ny-ny_im)/2]=pixelflux(xx,yy,params,cutoff=_EXTRA.cutoff) else $
  flux[(nx-nx_im)/2, (ny-ny_im)/2] = pixelflux(xx,yy,params,cutoff=0)

 ;test whether simply smoothing the function works, instead of convolving
 ; it does, suggesting the derivative may not be the biggest problem
 ;flux = smooth( flux, [4,4] )

 ;fft
im_fft = fft( flux, /double)
psf_fft = fft(psf, /double)
conv_fft = im_fft * psf_fft

flux = (double(fft(conv_fft,/double, /inverse)))*nx*ny
flux=rebin(flux,nx/ff,ny/ff)


return, flux[(nx-nx_im)/(2*ff):(nx-nx_im)/(2*ff)+nx_im/ff-1, $
             (ny-ny_im)/(2*ff):(ny-ny_im)/(2*ff) + ny_im/ff-1]

END


