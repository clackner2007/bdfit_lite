PRO oned_image, flux, radii, ivar, image, pos_ang, q, binwidth=binwidth, $
                xcen=xcen, ycen=ycen, numpixels=numpixels
;+
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
; May 21, 2010
; NAME: ONED_IMAGE
; calculates a 1-d profile of an image, using elliptical annuli
;
; INPUTS:  flux: output 1-D array of mean fluxes in each annulus
;          radii: radius of each flux bin
;          image
;          pos_ang = position angle (radians)
;          q = axis ratio (x/y) (big/small)
;          binwidth(optional) = width of each bin (pixels), default = 1
;          xcen, ycen (optional) center of ellipse, otherwise chosen
;          as maximum flux along each axis
;          numpixels = number of pixels in each flux bin
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;-
nx = (size(image, /dimensions))[0]
ny = (size(image, /dimensions))[1]

if not keyword_set(binwidth) then $
  binwidth = 1.0 ;max(sqrt((nx-xcen)*(ny-ycen))/30.,1.0)


if not keyword_set(xcen) then $
  xcen = (where(total(image,2) eq max(total(image,2))))[0]
if not keyword_set(ycen) then $
  ycen = (where(total(image,1) eq max(total(image,1))))[0]

;make array of distances to center along ellipese
dist_ellipse, ell, [nx, ny], xcen, ycen, q, pos_ang*180.0/!pi

annulus = ell le binwidth
j=where(annulus, np)
numpixels = [np]
flux = [total((image)*(annulus))/np]
radii = [binwidth-1.0]
ivar = [1.0/variance((image*annulus)[where(annulus ne 0)])]

num=1
while( 4*(num*binwidth)^2 lt nx*ny ) do begin
    annulus = (ell gt (num-1)*binwidth)*(ell le num*binwidth)
    j=where(annulus,np)
    if(np gt 0) then begin
        numpixels=[[numpixels],np]
        flux=[flux,[total(image*(annulus))/np]]
        
        ivar = [[ivar],1.0/variance((image*annulus)$
                                    [where(annulus ne 0)])]
    endif
    num+=1
    radii = [[radii],num*binwidth]
endwhile

;flux=[reverse(flux),flux[1:*]]
;radii=[reverse(radii),radii[1:*]]


END
