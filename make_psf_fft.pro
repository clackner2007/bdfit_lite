;;
; gets the fft of the psf at the proper image size
;;;

FUNCTION make_psf_fft, psfcrop, im

   xpsf = (size( psfcrop, /dimensions ))[0]
   ypsf = (size( psfcrop, /dimensions ))[1]
   nx_im = (size(im, /dimensions))[0]
   ny_im = (size(im, /dimensions))[1]
   nx = max([xpsf, nx_im])
   ny = max([ypsf, ny_im])

   psf = dblarr( nx, ny )

   psf[nx-xpsf/2, ny-ypsf/2] = psfcrop[0:(xpsf/2-1), 0:(ypsf/2-1)]
   psf[0, ny-ypsf/2] = psfcrop[xpsf/2:xpsf-1, 0:(ypsf/2-1)]
   psf[nx-xpsf/2, 0] = psfcrop[0:(xpsf/2-1), ypsf/2:ypsf-1]
   psf[0, 0] = psfcrop[xpsf/2:xpsf-1, ypsf/2:ypsf-1]

   return, fft(psf, /double)

END
