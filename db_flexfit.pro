;____________________________________
;
; Claire Lackner
;  Jan 11 2010
;
;  diskbulgefit.pro
;
;  fits an image with a bulge and a disk, the bulge with sersic
;  index=bulgeSersic and the disk with sersic index =1
; 
;  the disk parameters are returned in diskparams, and the bulge
;parameters in bulgeparams, the image, and the inverse variance are
;passed in via image and iv, the psf FWHM is given by psf, the
;chi-squared of the fit is returned by chisquare
;
;
; the fits have c=0 (bulge/disky isophotes) and the centers of the
;profiles are equal to each other
;____________________________________



PRO db_flexfit, params, image, psf, iv, chisquare, covar, $
                errors, status, degfree, skyVal=skyVal, $
                fix_params=fix_params, $
                free_sky=free_sky, $
                free_coords=free_coords, $
                negative=negative, fixq=fixq, $
                rescale=rescale, _EXTRA=_EXTRA
                

  xlen = (size(image, /dimensions))[0]
  ylen = (size(image, /dimensions))[1]

  if iv[0] eq -1 then begin
     good=where(image gt 1.0E-12 )
     iv = dblarr(xlen, ylen)
     iv[good] = 1.0/image[good]
;iv = dblarr(xlen,ylen)+1.0
  endif

  if( n_elements(psf) lt 2 ) then begin
     data = { image:image, $
              ivar:iv,$
              psf:dblarr(xlen,ylen) }
     data.psf = makepsf( xlen, ylen, psf )
  endif else data = {image:image, ivar:iv, psf:psf }
                                ;structure containing, image
                                ;and inverse variance and psf, 
                                ;make arrays for (x,y) values
  
  makegrid, xlen, ylen, x, y
  
;number of parameters, including sky
  nparams = n_elements(params)+1
;number of profiles
  nprof = n_elements(params)/8
  plist = indgen(nprof)*8
  phi_ind = 7
  x0_ind = 5
  bulge_ind = 2                 ;OR 10, figure this out
  skyInd = n_elements(params)

  if n_elements(fix_params) eq 0 then begin
     fix_params = intarr(nparams-1)
     fix_params[4+plist] = 1
  endif
  fixsky = 1
  if keyword_set(free_sky) then fixsky=0

  if (keyword_set(rescale) and (nprof gt 1)) then rescale=1 else rescale=0

  
  parinfo=replicate({value:0.D, fixed:0, $
                     limited:[0,0], limits:[0.D,0], $
                     mpside:2, tied:''}, nparams)
  
  if n_elements(skyVal) eq 0 then skyVal=0.0
  parinfo[*].value = [params, skyVal]
  parinfo[*].fixed = [fix_params, fixsky]

                                        ;define parinfo, which gives limits to
                                        ;parameters as well as starting values


  if (not keyword_set(free_coords)) then begin
     temp = x0_ind
     for prof=0, nprof-2 do begin
        match = x0_ind+8
        match2 = x0_ind+8+1
        parinfo[temp].tied = 'p['+string(match)+']'
        parinfo[temp+1].tied = 'p['+string(match2)+']'
        temp = match
     endfor
  endif


  if( not keyword_set(negative) ) then begin
     parinfo[plist].limited = [1,0]
     parinfo[plist].limits = [0.0, 0]
  endif
  parinfo[plist+1].limited = [1,0]
  parinfo[plist+1].limits = [0.0, sqrt(xlen*ylen*1.0)]
  if (rescale eq 1) then begin
     parinfo[plist[1:*]+1].limits=[1.0,1.e6] ;NOT an implemented upper limit, 
                                ;rescale assumes subsequent profiles
                                ;are larger than the first (BULGE
                                ;COMES FIRST!)
  endif

     parinfo[plist+2].limited = [1,1]
     parinfo[plist+2].limits = [0.1, 9.0] ; blanton uses 6.0
     reset = where(parinfo[plist+2].fixed ne 0, /null)
     if (reset ne !NULL) then $
        parinfo[(plist+2)[reset]].limits = $
        transpose([[parinfo[(plist+2)[reset]].value-1], $
                   [parinfo[(plist+2)[reset]].value+1]])
     
     parinfo[plist+3].limited = [1,1]
     parinfo[plist+3].limits = [0.05D,20.0D]
     reset = where(parinfo[plist+3].fixed ne 0, /null)
     if (reset ne !NULL) then $
        parinfo[(plist+3)[reset]].limits = $
        transpose([[parinfo[(plist+3)[reset]].value-1], $
                   [parinfo[(plist+3)[reset]].value+1]])

                                ;this is almost always fixed
     parinfo[plist+4].limited = [1,1]
     parinfo[plist+4].limits = [-1.9D, 1.9D]
     
     
                                ;coordinates should be on image
     parinfo[x0_ind+plist].limited = [1,1]
     parinfo[x0_ind+plist].limits = [0, xlen+0.5]
     parinfo[x0_ind+plist+1].limited = [1,1]
     parinfo[x0_ind+plist+1].limits = [0, ylen+0.5]


  start_params = parinfo[*].value

                                        ;normalize psf flux to one
  ;data.psf /= total(data.psf)

                                        ;normalize maximum pixel value
                                        ;to maximum pixel flux
  tot = max([total(data.image),10])
  ;; if( nparams lt 16) then $
  ;;    parinfo[0].value = tot/totalsersicflux(start_params) $
  ;; else begin
  ;;    rb=start_params[9]
  ;;    if rescale ne 0 then rb *= start_params[1]
  ;;    d = max([totalsersicflux([start_params[0:4],0.0,0.0,0.0]),1.0e-2])
  ;;    b = max([totalsersicflux([start_params[8], $
  ;;                              rb, $
  ;;                              start_params[10:12],0.0,0.0,0.0]),1.0e-2])
     
  ;;    parinfo[8].value = frac*tot / b
  ;;    parinfo[0].value = tot*(1.0-frac)/d
  ;; endelse
  ;; start_params = parinfo[*].value
                                        ;DON'T USE errors, from DATA.ERR,
                                        ;these are inverse variances, not
                                        ;errors, rely on weights keyword for
                                        ;errors!!!
  temp = start_params


  again=1
  fitstat=1
  times = 0
  reset = 0
  seed=systime(/seconds)

  while( again ne 0 and times lt 5 ) do begin
    
     start_params = temp
     parinfo[*].value = start_params
     again = 0
     fitstat = 1
     print, start_params
     params=mpfit2dfun('pixelfluxpsf', x, y, data.image, $
                       data.ivar, $
                       start_params, parinfo=parinfo, $
                       functargs={$;psfFFT:data.psf,  $
                                  psfImage:data.psf, $
                                  cutoff:1, rescale:rescale}, $
                       perror=errors, covar=covar, weights=data.ivar,$
                       dof=dof, bestnorm=bn, $ ;/quiet, $
                       xtol=1.0e-10, gtol=1.0e-8, ftol=1.0e-6, $
                       status=fitstat, maxiter=200, npegged=npg) ;, $
                                ;iterproc='iter_plot', ITERARGS={psf:data.psf}, yfit=fit)

     
    
    if( nparams ge 16 and times lt 4 ) then begin
        ;don't allow one set of parameters to be zero
        if (params[1] - parinfo[1].limits[0]) lt 1.e-8 then begin
            params[1] = start_params[1]*(1.-0.2*times)
            again = 1                
         endif
        if parinfo[1].limited[1] eq 1 then begin
           if (parinfo[1].limits[1] - params[1]) lt 1.e-3 then begin
              if rescale ne 0 then params[9] = params[9]*params[1]
              params[1] = (0.9-0.1*times)*start_params[1]
              again = 1
              if rescale ne 0 then params[9] = params[9]/params[1]
           endif
        endif
        if ((params[3] - parinfo[3].limits[0]) lt 1.e-8) then begin
            params[3] = (start_params[3]-parinfo[3].limits[0])*$
              (1.-0.2*times) + parinfo[3].limits[0]*1.1
            again=1
        endif
        if (params[0] - parinfo[0].limits[0]) lt 1.0e-8 then begin
            d=max([totalsersicflux([1.0,params[1:4],0.0,0.0,0.0]),1.0e-2])
            params[0] = (1.0-frac)*tot/(d)
            again = 1
        endif
        ;don't allow for no bulge
        if (params[9] - parinfo[9].limits[0]) lt 1.e-8 then begin
            params[9] = start_params[9]*(1.-0.2*times)
            again=1
        endif
        ;; if (rescale ne 0) and $
        ;;    (parinfo[9].limits[1] - params[9] lt 1.0e-5) then begin
        ;;    params[9] = (1.9-0.1*times)
        ;;    again=1
        ;; endif
        if(params[11] - parinfo[11].limits[0]) lt 1.e-8 then begin
            params[11] =  (start_params[11]-parinfo[11].limits[0])*$
              (1.-0.2*times) + parinfo[11].limits[0]*1.1
            again=1
        endif
        if (params[8] - parinfo[8].limits[0]) lt 1.0e-8 then begin
           rb = params[9]
           if rescale ne 0 then rb *= params[1]
            b=max([totalsersicflux([1.0,rb, $
                                    params[10:12],0.0,0.0,0.0]),1.0e-2])
            params[8] = 0.5*tot/b
            again = 1
        endif
 
    endif 
    ;make q <=1 for disk and bulge
    for prof=0, nprof-1 do begin
       curr = prof*8
       if( (params[curr+3] gt 1.0) and abs(params[curr]) gt 1.0e-18 $
           and (params[curr+1] gt 1.0e-12) ) then begin
          params[curr+3] = 1.0/params[curr+3]
          params[curr+phi_ind] += !pi*0.5
          params[curr+1] /= params[curr+3]
          if (rescale ne 0 ) then begin
             if (nprof gt 1) and (curr eq 0) then begin
                for pp=1, nprof-1 do $
                   params[curr+8*pp+1] *= params[curr+3]
             endif 
             if (curr gt 0) and (params[curr+1] lt 1.0) then begin
                again = 1
                params[curr+1] = 1.2
             endif
          endif
          if (again eq 1) then $
             params[curr+3] = max([params[curr+3], 1.1*parinfo[curr+3].limits[0]])
          if (params[curr+1] ge parinfo[curr+1].limits[1]) then begin
             params[curr+1] = start_params[curr+1]*(1-0.2*times)
             again = 1
          endif
       endif
       params[curr+phi_ind] = params[curr+phi_ind] mod !pi
    endfor
    
    ;if the fit didn't work, try again
    if( fitstat eq 5 ) then begin
        again=1
        params[0] *= 0.9
        if( nparams ge 16) then begin
            params[8] *= 1.1
        endif
     endif
    
    temp = params
    times += 1
 endwhile

  if(times ge 5 and again eq 1) then begin
     print, 'reached max number of tries before fit was found'
     print, 'status=', fitstat
     print, 'times=',times
     print, 'params=',params
     params[*] = -1.0
    fitstat = 0
 endif


;put the rotation in units of pi
  
  params[plist+phi_ind] = params[plist+phi_ind] mod !pi

  if keyword_set(rescale) then begin
     errors[plist[1:*]+1] = sqrt(params[plist[1:*]+1]^2*errors[plist[1:*]+1]^2 + $
                                 params[1]^2*errors[1]^2 + $
                                 2*params[plist[1:*]+1]*params[1]*covar[1,plist[1:*]+1])
     covar[*,plist[1:*]+1] *= params[1] ;THIS IS WRONG!!!
     covar[plist[1:*]+1,*] *= params[1] ;still WRONG
     params[plist[1:*]+1] *= params[1]
  endif
;don't count null pixels in the final DOF to get the reduced chi squared
  junk=where(data.ivar lt 1.0e-12, nullpix)
  chisquare=bn/(dof - nullpix)
  status=fitstat
  degfree = dof-nullpix
  skyVal = params[skyInd]
  params = params[0:skyInd-1]

END
