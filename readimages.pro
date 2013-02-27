;____________________________________
;
; Claire Lackner
;  Sept 7, 2009
;
;  readimages.pro
;
;  IDL
; reads the psf, noise and image from sdss dr4 lowz catalog
; returns a structure with image, ivar, and psf
;____________________________________

;function takes object name as argument and desired band ugriz=(01234)

FUNCTION readimages, folder, name, band, status=status

bands=['u','g','r','i','z']


index = where( bands[*] eq band )

im_hdu = index[0]
ivar_hdu = index[0]*2+1
psf_hdu = index[0]

;read images from dr6 stuff, but keep psf from here
path='/peyton/scr/depot1/clackner/vagc2/lowz/images/'+folder+'/'
imfile=path+'lowz-atlas-'+name+'.fits'
psfile=path+'psf-'+name+'.fits'
path6='/peyton/scr/depot1/clackner/vagc2/lowz_dr6/images/'+folder+'/'
imfile6=path6+name+'-child-fits'
errfile6=path6+name+'-parent-fits'

im=mrdfits(imfile6, im_hdu, status=imstat, /silent)
ivar=mrdfits(errfile6, ivar_hdu, status=varstat, /silent)

;get the psf from the full width half-max in the data
;psf=mrdfits(psfile, psf_hdu, status=psfstat, /silent)
;but keep a place for the psf in data structure
psf = im*0

if( imstat lt 0 or varstat lt 0 ) then begin
    print, "can't open file: "+name
    status = -1
endif else status = 1

;err=sqrt(1./ivar)
data={image:im, ivar:ivar, psf:psf }

return,data

end
