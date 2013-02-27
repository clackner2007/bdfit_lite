PRO onedfit, diskparams, bulgeparams, flux, psf_fwhm, iv, $
             chisquare, covar, errors, status, degfree, bulgeSersic=bulgeSersic, $
             nodisk=nodisk, freebulge=freebulge, $
             negative=negative, q=q, pa=pa
;+
;;;;;;;;;
; Claire Lackner
; May 21 2010
;
; onedfit.pro
; program for fitting one-d profiles of images
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;-


x = dindgen(n_elements(flux))

if( not keyword_set(nodisk) ) then begin
    nparams = 14
    phi_ind = 5
    x0_ind = 12
    bulge_ind = 8
endif else begin
    nparams = 8
    phi_ind = 7
    x0_ind = 5
    bulge_ind = 2
endelse
parinfo=replicate({value:0.D, fixed:0, $
                   limited:[0,0], limits:[0.D,0]}, nparams)

if(keyword_set(nodisk)) then parinfo[*].value = [1.0,10.,bulgeSersic, $
                                                 1.0,0.0,0.0,0.0,0.0] $
else parinfo[*].value = [1.0, 15.1, 1.0, q, $
                         0.0, pa, 1.0, 15.1, bulgeSersic, q, 0.0, pa, $
                         n_elements(flux),n_elements(flux) ]

if( not keyword_set(negative) ) then begin
    parinfo[0].limited = [1,0]
    parinfo[0].limits = [0.0, 0]
endif
parinfo[1].limited = [1,0]
parinfo[1].limits = [0.1,0.0]  ;smaller than the seeing, but not so
                               ;small to get ewrid limits
parinfo[2].limited = [1,1]
parinfo[2].limits = [0.05, 6.0]       ; blanton uses 6.0
parinfo[2].fixed = 1
parinfo[3].limited = [1,0]
parinfo[3].limits = [0.01D,0]
parinfo[4].limited = [1,1]
parinfo[4].limits = [-1.9D, 1.9D]
parinfo[4].fixed = 1
if( nparams eq 14 ) then begin
    if not keyword_set(negative) then begin
        parinfo[6].limited = [1,0]
        parinfo[6].limits = [0.0, 0]
    endif
    parinfo[7].limited = [1,0]
    parinfo[7].limits = [0.1,0.0]
    parinfo[8].limited = [1,1]
    parinfo[8].limits = [0.05, 6.0]     ; blanton uses 6.0
    parinfo[8].fixed = 1
    parinfo[9].limited = [1,0]
    parinfo[9].limits = [0.01D,0]
    parinfo[10].limited = [1,1]
    parinfo[10].limits = [-1.9D, 1.9D]
    parinfo[10].fixed = 1
endif
parinfo[x0_ind].limited = [1,1]
parinfo[x0_ind].limits = [0.0, 2*n_elements(flux)]
parinfo[x0_ind+1].limited = [1,1]
parinfo[x0_ind+1].limits = [0.0, 2*n_elements(flux)]

if (keyword_set(freebulge)) then parinfo[bulge_ind].fixed = 0             

start_params = parinfo[*].value

norm= max(flux) / $
  max(pixelfluxpsf1d(x,start_params, _EXTRA={psfImage:psf_fwhm,$
                                             cutoff:1,q:q, pa:pa}))
if keyword_set(nodisk) then begin
    parinfo[0].value = norm
    start_params[0] = norm
endif else begin
    parinfo[0].value = norm
    parinfo[6].value = norm
    start_params[0] = norm
    start_params[6] = norm
endelse
       

params=mpfitfun('pixelfluxpsf1d', x,flux,ivar,$
                start_params, parinfo=parinfo, $
                functargs={cutoff:1,psfImage:psf_fwhm,q:q, pa:pa}, $
                dof=dof, bestnorm=bn, $;/quiet, $
                xtol=1.0e-8, gtol=1.0e-8, ftol=1.0e-8, $
                status=fitstat, maxiter=75 )

print, params
diskparams = params

END

