;____________________________________
;
; Claire Lackner
;  Jan 11 2010
;
; coverage.pro
;
;  creates a circle in the center of the image which has a size given
;  by scale, relative to the size of the image
;____________________________________


FUNCTION coverage, ivar, scaling, xcenter, ycenter, radius

nx = (size( ivar, /dimensions))[0]
ny = (size( ivar, /dimensions))[1]

w = dblarr( nx, ny )



;weights zero outside central region
w = make_circ( w, scaling*radius, cent=[xcenter, ycenter] ) * ivar

return, w

END
