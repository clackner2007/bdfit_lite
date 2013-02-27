;____________________________________
;
; Claire Lackner
;  Jan 11 2010
;
;  fitfakes.pro
;
;  like sersicfit.pro
;  fits fake images of galaxies created by makegal.pro
; 
; the output is a comparison of the input parameters
; to the fitted parameters
; 
;____________________________________



PRO fitfakes, profile, vals, fit_range, image

;data is (image, noise, psf) made public so function can access psf
;common datastruct, data

psf = abs(randomn(seed)*0.5 + 2.6)
xlen = 0
ylen = 0
inputs = 0
data = { image:makefake(profile, xlen, ylen, inputs, psf=psf), $
         ivar:dblarr( xlen, ylen ) + 1.0D, $
         psf:dblarr(xlen,ylen) }

                                        ;structure containing, image
                                        ;and inverse variance and psf, 

                                        ;crop unnecessary noise on outside of image
newdata = cropimage( data.image, data.ivar, max([xlen, ylen]), $
                     x0=xcrop, y0=ycrop )

                                        ;combine new and old data
data = create_struct( data, newdata )

imsize = size(data.imagecrop, /dimensions)

                                        ;make arrays for (x,y) values
xrow = dindgen( imsize[0] )
ycol = dindgen( imsize[1] )

x = xrow#(ycol*0+1)
y = (xrow*0+1)#ycol

                                        ;center model at highest value pixel
indmax = where( data.imagecrop[*] eq max(data.imagecrop) )
y0 = indmax[0] / imsize[0]
x0 = indmax[0] mod imsize[0]

                                        ;make the psf image, if it isn't there yet
if( total(data.psfcrop) lt 0.04 ) then $
  data.psfcrop = makepsf( imsize[0], imsize[1], $
                          psf )

parinfo=replicate({value:0.D, fixed:0, $
                   limited:[0,0], limits:[0.D,0]}, 8)

                                        ;define parinfo, which gives limits to
                                        ;parameters as well as starting values

;vals contains the original input parameters, as well as the fit at
;each size
if( profile eq 1 ) then n_prof = 1 else n_prof = 2
vals = dblarr( n_prof + n_elements(fit_range), 10 )
vals[0:n_prof-1,0:7] = transpose(inputs)
vals[0:n_prof-1,8] = psf + dblarr(n_prof)

;do the fitting from scratch for each value in fit_range
for s=0,n_elements(fit_range)-1 do begin

    parinfo[*].value = [1.0, 21.0, 3.1, $
                        0.9, $
                        -0.1, $
                        x0, y0, 0.78 ]

    parinfo[0].limited = [1,0]
    parinfo[0].limits = [1D-8, 0]
    parinfo[2].limited = [1,1]
    parinfo[2].limits = [1D-8, 14.0]    ; blanton uses 6.0
    parinfo[1].limited = [1,0]
    parinfo[1].limits = [1D-8,0]
    parinfo[3].limited = [1,0]
    parinfo[3].limits = [0.0D,1.0000001D]
    parinfo[4].limited = [1,1]
    parinfo[4].limits = [-1.9D, 1.9D]
    parinfo[5].limited = [1,1]
    parinfo[5].limits = [0.0, imsize[0]+0.05]
    parinfo[6].limited = [1,1]
    parinfo[6].limits = [0.0, imsize[1]+0.05]

    start_params = parinfo[*].value


                                        ;normalize psf flux to one
    data.psfcrop /= total(data.psfcrop)

                                        ;normalize maximum pixel value
                                        ;to maximum pixel flux
    norm = max(data.imagecrop) / $
      max(pixelfluxpsf(x, y, start_params, _EXTRA={psfImage:data.psfcrop}))

    parinfo[0].value = norm
    start_params[0] = norm
                                        ;DON'T USE errors, from DATA.ERR,
                                        ;these are inverse variances, not
                                        ;errors, rely on weights keyword for
                                        ;errors!!!

    temp = start_params
    again=1
    times = 0
    fitstat = 0
    while( again ne 0 and times lt 10 ) do begin
        
        start_params = temp
        again = 0
        
        params=mpfit2dfun('pixelfluxpsf', x, y, data.imagecrop, $
                          data.ivarcrop, $
                          start_params, parinfo=parinfo, $
                          weights = coverage(data.ivarcrop, fit_range[s]), $
                          perror=errors, covar=covar, $
                          dof=dof, bestnorm=bn, /quiet, status=fitstat, $
                          functargs={psfImage:data.psfcrop})
        
        
        if( params[3] gt 1.0 ) then begin
            again = 1
            params[3] = 1.0/params[3]
            params[7] -= !PI/2.0D
        endif

        temp = params
        times += 1
    endwhile

    if(fitstat eq 5 or times gt 9) then params[*] = 0.0


    params[5] += xcrop
    params[6] += ycrop

    vals[s+n_prof, 0:7] = params
    vals[s+n_prof,8] = fit_range[s]
    vals[s+n_prof,9] = bn/dof
  
endfor

image = data.image

END


