;____________________________________
;
; Claire Lackner
;  Jan 11 2010
;
;  fitImage.pro
;
;  fits an image at certain values of the measured half-lgiht radius
;  and put those values in vals, including a fit to the entire image
; 
;____________________________________



PRO fitimage, vals, fit_range, image, psf, iv


xlen = (size(image, /dimensions))[0]
ylen = (size(image, /dimensions))[1]


if( iv[0] lt 0.5 ) then begin
; good=where(image gt 0.0 )
; iv = dblarr(xlen, ylen)
; iv[good] = 1.0/image[good]
    iv = dblarr(xlen, ylen)+1.0D
endif

data = { image:image, $
         ivar:iv,$
         psf:dblarr(xlen,ylen) }
                                        ;structure containing, image
                                        ;and inverse variance and psf, 
                                        ;make arrays for (x,y) values
xrow = dindgen( xlen )
ycol = dindgen( ylen )

x = xrow#(ycol*0+1)
y = (xrow*0+1)#ycol

                                        ;center model at highest value pixel
indmax = where( data.image[*] eq max(data.image) )
y0 = indmax[0] / xlen
x0 = indmax[0] mod xlen

                                        ;make the psf image, if it isn't there yet
if( total(data.psf) lt 0.04 ) then $
  data.psf = makepsf( xlen, ylen, psf )

parinfo=replicate({value:0.D, fixed:0, $
                   limited:[0,0], limits:[0.D,0]}, 8)

                                        ;define parinfo, which gives limits to
                                        ;parameters as well as starting values

;the fit values go in vals
vals = dblarr( n_elements(fit_range)+1, 11 )
vals[1:*,8] = fit_range

parinfo[*].value = [1.0, 10.1, $
                    3.8, $
                    0.85, $
                    0.0, $
                    x0, y0, 0.2 ]

parinfo[0].limited = [1,0]
parinfo[0].limits = [1D-8, 0]
parinfo[1].limited = [1,0]
parinfo[1].limits = [1D-3,0.0]
parinfo[2].limited = [1,1]
parinfo[2].limits = [1D-3, 8.0]         ; blanton uses 6.0
parinfo[3].limited = [1,0]
parinfo[3].limits = [0.001D,1.0000001D]
parinfo[4].limited = [1,1]
parinfo[4].limits = [-1.9D, 1.9D]
parinfo[4].fixed = 1
parinfo[5].limited = [1,1]
parinfo[5].limits = [0.0, xlen+0.05]
parinfo[6].limited = [1,1]
parinfo[6].limits = [0.0, ylen+0.05]

start_params = parinfo[*].value

                                        ;normalize psf flux to one
data.psf /= total(data.psf)

                                        ;normalize maximum pixel value
                                        ;to maximum pixel flux
norm = max(data.image) / $
  max(pixelfluxpsf(x, y, start_params, $
                   _EXTRA={psfImage:data.psf}))

parinfo[0].value = norm
start_params[0] = norm
                                        ;DON'T USE errors, from DATA.ERR,
                                        ;these are inverse variances, not
                                        ;errors, rely on weights keyword for
                                        ;errors!!!
temp = start_params

;do the fitting from scratch for each value in fit_range
for s=0,n_elements(fit_range) do begin
    again=1
    fitstat=1
    times = 0

;    if( s ne 0 ) then temp[1] = vals[0,1]*0.5
;    if temp[2] gt 7.99 then temp[2] = 4.1
    
    while( again ne 0 and times lt 4 ) do begin
        
        start_params = temp
        parinfo[*].value = start_params
        again = 0
        fitstat = 1

        if( s eq 0 ) then begin
            params=mpfit2dfun('pixelfluxpsf', x, y, data.image, $
                              1.0/sqrt(data.ivar), $
                              start_params, parinfo=parinfo, $
                              functargs={psfImage:data.psf}, $
                              perror=errors, covar=covar, weights=data.ivar,$
                              dof=dof, bestnorm=bn, /quiet, $
                              status=fitstat, maxiter=100 )
        endif else begin
            params=mpfit2dfun('pixelfluxpsf', x, y, data.image, $
                              data.ivar, $
                              start_params, parinfo=parinfo, $
                              functargs={psfImage:data.psf}, $
                              weights = coverage(data.ivar, $
                                                 fit_range[s-1], $
                                                 vals[0,5], $ ;xcenter
                                                 vals[0,6],$  ;ycenter
                                                 vals[0,1]),$ ;halflightradius
                              perror=errors, covar=covar, $
                              dof=dof, bestnorm=bn, /quiet, $
                              status=fitstat, maxiter=100 )
        endelse
        
        if( params[3] gt 1.0 ) then begin
            again = 1
            params[3] = 1.0/params[3]
            params[7] -= !PI/2.0D
        endif

        temp = params
        times += 1
    endwhile

    if(fitstat eq 5 or times gt 4) then begin
        params[*] = 0.0
        print, 'reached max number of tries before fit was found'
        print, 'status=', fitstat
        print, 'times=',times
    endif
        
    vals[s, 0:7] = params
    vals[s,9] = bn
    vals[s,10] = dof
    
endfor


END


