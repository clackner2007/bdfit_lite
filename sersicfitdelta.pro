;____________________________________
;
; Claire Lackner
;  Sept 7, 2009
;
;  sersicfitdelta.pro
;
;  IDL
; gives the bestfit sersic model for a galaxy from sdss low-z catalog,
; assuming a delta function psf, the data is given in the files  from
; readgal.pro 
; 
; the fitting is done using mpfitfun.pro (citation) the sersic
; function is given by pixelflux.pro, which gives the flux integrated
; over a pixel
;   
;____________________________________


PRO sersicfitdelta

common datastruct, data

gals = lowz()

for i=0L, n_elements(gals.name)-1L do begin
    ;structure containing, image, err, psf
     data = readimages( gals.name[i], 'r' )

    ;ignore the psf in this program
    imsize = size(data.image, /dimensions)
    
    xrow = dindgen( imsize[0] )
    ycol = dindgen( imsize[1] )
    
    x = xrow#(ycol*0+1)
    y = (xrow*0+1)#ycol

    indmax = where( data.image[*] eq max(data.image) )
    y0 = indmax[0] / imsize[0]
    x0 = indmax[0] mod imsize[0]

    parinfo=replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 8)


    parinfo[*].value = [1.0, gals.sersic_n[i], $
                        gals.reff[i], 0.5, -0.1, x0, y0, 0.78 ]

    parinfo[0].limited = [1,0]
    parinfo[0].limits = [0.0, 0]
    parinfo[1].limited = [1,1]
    parinfo[1].limits = [0.0, 14.0]
    parinfo[2].limited = [1,0]
    parinfo[2].limits = [0.0,0]
    parinfo[3].limited = [1,1]
    parinfo[3].limits = [0.0,1.00001]
    parinfo[5].limited = [1,1]
    parinfo[5].limits = [0.0, imsize[0]+0.05]
    parinfo[6].limited = [1,1]
    parinfo[6].limits = [0.0, imsize[1]+0.05]

    start_params = parinfo[*].value


    norm = max(data.image) / max(pixelflux(x, y, start_params))

    parinfo[0].value = norm
    start_params[0] = norm
    params=mpfit2dfun( 'pixelflux', x, y, data.image, data.err, $
                       start_params, parinfo=parinfo, $ 
                       perror=errors, covar=covar, dof=dof, bestnorm=bn )

;    params = start_params
    print, gals.name[i], params, errors, dof, bn
    atv, pixelflux( x, y, params) - data.image

endfor
END

