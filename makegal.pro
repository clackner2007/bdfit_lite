;____________________________________
;
; Claire Lackner
;  Jan 11, 2009
;
;  makegal.pro
;
;  creates a centered image of galaxy with given
;  sersic profiles, all parameter measurements should be
;  given in pixels
;
;   These images can be PSF convolved
;
;   the size of the final image is either given by x, y, or it's
;   set automatically by the size of the galaxy
;
;   image = makegal(n_prof, params[8, n_prof], [, psf=psf_fwhm,
;   size=[x,y]])
;
;INPUTS:
;  params = n_prof x 8 array of Sersic parameters, with lengths given
;  in pixels
;  profile = This is the following for different galaxies
;      Elliptical:  profile = 0
;      Disk+Bulge:  profile = 1
;      Disk+disk:   profile = 2
;OPTIONAL KEYWORDS:
;  psf = FWHM of psf if a convolution is desired, if unset this step
;  isn't done
;  size = dimensions of final 2d array, if not set this is
;  automatically chosen
;
;OUTPUT:
; 2-d array of centered galaxy image
;____________________________________

;***************************************
;function that returns the flux of the galaxy when it is a sum of
;Sersic profiles and disk profile, the size of the image is given by
;nx and ny, params is an 8x2 or 8x1 array, with the disk profile
;always coming first
 FUNCTION fluxsum, x, y, profile, params

 im = dblarr( (size(x, /dimensions))[0], (size(x, /dimensions))[1] )
 ; only a single sersic profile
 if( profile eq 1 ) then im = sersic( x,y, params )
 if( profile eq 2 ) then im = disk(x,y,params[*,0]) + sersic(x,y,params[*,1])
 if( profile eq 3 ) then im = disk(x,y,params[*,0]) + disk(x,y,params[*,1])

 return, im

 END

;***************************************
FUNCTION makegal, profile, params, psf=psf, imdimen=imdimen

if( n_params() ne 2 ) then message, $
  'usage: makegal(type of profile, params[n_profiles, 8], [psf=fwhm, imdimen=[x,y]])'

if not keyword_set(imdimen) then $
  imdimen = [ max( [120, fix(3 * max( params[2, *]  ))] ), $
              max( [120, fix(3 * max( params[2, *] ))] ) ]

image = dblarr( imdimen[0], imdimen[1] )
nx = imdimen[0]
ny = imdimen[1]
xx = dindgen( nx )
yy = dindgen( ny )
px = xx#(yy*0+1)
py = (xx*0+1)#yy
x0 = px[0]
y0 = py[0]


flux0 = fluxsum( px, py, profile, params )

i=1

THRESH = 1.0e-2
bad = [1]
while i lt 6 and bad[0] ne -1 do begin
    row = (0.5+dindgen(nx*2^i))/(2^i) + x0 - 0.5
    col = (0.5+dindgen(ny*2^i))/(2^i) + y0 - 0.5

    x = row#(col*0+1.)
    y = (row*0+1)#col
    
    flux_large = fluxsum( x, y, profile, params )
    fluxbin = rebin( flux_large, nx, ny )
    
                                        ;put the binned flux into the old flux 
    fluxnew = flux0
    fluxnew[x0,y0] = fluxbin
    
                                        ;take the difference and selecte the pixels that didn't converge
    notzeros = where( flux0 gt 1.0e-15, count )

    difflux = abs( flux0 - fluxnew )
    if( count gt 0 ) then difflux[notzeros] = difflux[notzeros]/flux0[notzeros]
    bad = where(difflux gt THRESH)
    
                                        ;get out of while loop if there are no elements under the 
                                        ;threshold
    if bad[0] eq -1 then continue
    
                                        ;set the new parameters for the area that need more subdividing
    x0 = min( px[bad] )
    y0 = min( py[bad] )
    nx = max( px[bad] ) - x0 + 1
    ny = max( py[bad] ) - y0 + 1

    flux0 = fluxnew

 i++
endwhile

;if there is a psf, do the convolution
if keyword_set(psf) then begin
    nx = imdimen[0]
    ny = imdimen[1]
                                  
    psfImage=makepsf(nx, ny, psf )
                                        ;the psf needs to be divided
                                        ;into quadrants and rearranged
                                        ;to get the correct image
    psfScramble=dblarr(nx, ny)
    
; psf, with the 'negtive' frequencies placed in the opposite edges
    psfScramble[nx-nx/2, ny-ny/2] = psfImage[0:(nx/2-1), 0:(ny/2-1)]
    psfScramble[0, ny-ny/2] = psfImage[nx/2:nx-1, 0:(ny/2-1)]
    psfScramble[nx-nx/2, 0] = psfImage[0:(nx/2-1), ny/2:ny-1]
    psfScramble[0, 0] = psfImage[nx/2:nx-1, ny/2:ny-1]

    ;fft
    im_fft = fft( flux0, /double )
    psf_fft = fft( psfScramble, /double )
    conv_fft = im_fft * psf_fft
    flux0 = (double(fft(conv_fft, /double, /inverse)))*imdimen[0]*imdimen[1]
endif

return, flux0

END
