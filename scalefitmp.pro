FUNCTION scalefitmp, image, psf, ivar, params
imsize = size(image, /dimensions)
makegrid, imsize[0], imsize[1], x,y

if( n_elements(params) eq 8 ) then begin
    parinfo=replicate({value:0.D, fixed:1, limited:[0,0], $
                       limits:[0.0D,0]},8)
    
    parinfo[*].value = params
    parinfo[0].fixed=0
    parinfo[0].limited=[1,0]
    parinfo[0].limits=[1.0e-8,0]
    
endif else begin
    parinfo=replicate({value:0.D, fixed:1, limited:[0,0], $
                       limits:[0.0D,0]},14)
    
    parinfo[*].value = params
    parinfo[0].fixed=0
    parinfo[6].fixed=0
    parinfo[0].limited=[1,0]
    parinfo[0].limits=[1.0e-9,0]
    parinfo[6].limited=[1,0]
    parinfo[6].limits=[1.0e-9,0]
endelse
start_params=params

newparams = mpfit2dfun('pixelfluxpsf', x,y,image, $
                       0.0, start_params, parinfo=parinfo,$
                       functargs={psfImage:psf, cutoff:1},$
                       perror=errors,covar=covar,weights=ivar,$
                       dof=dof,bestnorm=bn,/quiet,$
                       status=fitstat)
print, 'chisquare, mp ', bn, dof
print, 'mp ', newparams[0], newparams[6]
return, newparams

END
