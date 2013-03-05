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


; diskparams=output for disk parameters (in sersic form)
;bulgeparams=output for bugle parameters (
PRO diskbulgefit, diskparams, bulgeparams, image, psf, iv, $
                  chisquare, covar, errors,  status, degfree, $
                  skyVal, free_sky=free_sky, $
                  bulgeSersic=bulgeSersic, $
                  nodisk=nodisk, freebulge=freebulge, $
                  freedisk=freedisk, $
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

                                        ;center model at highest
                                        ;summed pixel in each
                                        ;direction (as in cropimage.pro)
;y0 = (where(total(data.image,1) eq max(total(data.image,1))))[0]
;x0 = (where(total(data.image,2) eq max(total(data.image,2))))[0]
ys=reverse(sort(smooth(total(data.image,1),min([4, ylen-1]))))
xs=reverse(sort(smooth(total(data.image,2),min([4, xlen-1]))))
y0=ys[0]
x0=xs[0]


if( not keyword_set(nodisk) ) then begin
    nparams = 17
    phi_ind = 7
    x0_ind = 5
    bulge_ind = 10
    skyInd = 16
endif else begin
   rescale=0
    nparams = 9
    phi_ind = 7
    x0_ind = 5
    bulge_ind = 2
    skyInd = 8
endelse

parinfo=replicate({value:0.D, fixed:0, $
                   limited:[0,0], limits:[0.D,0], $
                   mpside:2, tied:''}, nparams)

                                        ;define parinfo, which gives limits to
                                        ;parameters as well as starting values

if(keyword_set(nodisk)) then parinfo[*].value = [1.0,10.,bulgeSersic, $
                                                0.7,0.0,x0,y0,0.2, 0.0] $
else begin
    parinfo[*].value = [1.0, 30.1, 1.0, 0.85, $
                        0.0, x0, y0, 0.2, $
                        1.0, 8.1, bulgeSersic, 0.7, 0.0, $
                        x0, y0, 0.2, 0.0]
    
    if (not keyword_set(free_coords)) then begin
        parinfo[x0_ind].tied = 'p[13]'
        parinfo[x0_ind+1].tied = 'p[14]'
     endif
    if (keyword_set(rescale)) then begin
       rescale=1
       parinfo[9].value = 8.1/30.1
    endif else begin
       rescale=0
    endelse
    ;tie the bulge ellipticity to the disk ellipticity (and PA)
    ;parinfo[11].tied = 'p[3]'
    ;parinfo[15].tied = 'p[7]'
    

endelse 
frac=0.5

if( keyword_set(_EXTRA) ) then begin

    if tag_exist(_EXTRA, 'start') then parinfo[*].value=_EXTRA.start $
    else begin
        
        if(tag_exist(_EXTRA,"Reff")) then begin
            parinfo[1].value = _EXTRA.Reff
            if(not keyword_set(nodisk)) then begin
               if rescale eq 0 then $
                  parinfo[9].value = _EXTRA.Reff $
               else parinfo[9].value = 0.8
                frac=0.5
                if tag_exist(_EXTRA,"fracdev") then begin
                    frac=(0.10)+_EXTRA.fracdev*0.8
                    if(_EXTRA.fracdev lt 0.5) then $
                      parinfo[9].value *= 0.5 $
                    else $
                      parinfo[1].value *= 1.0 ;0.5 + _EXTRA.fracdev ;-1.+2.*_EXTRA.fracdev
                endif

            endif 
        endif 
        if(tag_exist(_EXTRA,"q")) then begin
            parinfo[3].value = max([0.11,_EXTRA.q])
            if(not keyword_set(nodisk)) then parinfo[11].value = max([0.11,_EXTRA.q])
        endif
        if(tag_exist(_EXTRA,"phi")) then begin
            parinfo[phi_ind].value = _EXTRA.phi
            if(not keyword_set(nodisk)) then parinfo[15].value = _EXTRA.phi
        endif
    endelse
endif 

if( not keyword_set(negative) ) then begin
    parinfo[0].limited = [1,0]
    parinfo[0].limits = [0.0, 0]
endif
parinfo[1].limited = [1,0]
parinfo[1].limits = [0.0, 2.0*sqrt(xlen*ylen*1.0)] ;smaller than the seeing, but not so
                               ;small to get ewrid limits
parinfo[2].limited = [1,1]
parinfo[2].limits = [0.1, 9.0]       ; blanton uses 6.0
parinfo[2].fixed = 1
if keyword_set(fixq) then begin
    parinfo[3].fixed=1
    if(tag_exist(_EXTRA,"q")) then $
      parinfo[3].value = _EXTRA.q $
    else $
      parinfo[3].value = 1.0
endif else begin
parinfo[3].limited = [1,1]
parinfo[3].limits = [0.05D,20.0D]
endelse 
parinfo[4].limited = [1,1]
parinfo[4].limits = [-1.9D, 1.9D]
parinfo[4].fixed = 1
if( nparams ge 16 ) then begin
    if not keyword_set(negative) then begin
        parinfo[8].limited = [1,0]
        parinfo[8].limits = [0.0, 0]
     endif
    parinfo[9].limited = [1,0]
    parinfo[9].limits = [0.0,1.01]
    if keyword_set(rescale) then begin 
       parinfo[9].limited=[1,1]
       parinfo[9].limits=[0.0,1.0]
    endif
    parinfo[10].limited = [1,1]
    parinfo[10].limits = [0.1, 9.0]     ; blanton uses 6.0
    parinfo[10].fixed = 1
    parinfo[11].limited = [1,1]
    parinfo[11].limits = [0.05D,20.0D]
    parinfo[12].limited = [1,1]
    parinfo[12].limits = [-1.9D, 1.9D]
    parinfo[12].fixed = 1
    parinfo[x0_ind+8].limited = [1,1]
    parinfo[x0_ind+8].limits = [0.0, xlen+0.05]
    parinfo[x0_ind+9].limited = [1,1]
    parinfo[x0_ind+9].limits = [0.0, ylen+0.05]
endif
parinfo[x0_ind].limited = [1,1]
parinfo[x0_ind].limits = [0.0, xlen+0.05]
parinfo[x0_ind+1].limited = [1,1]
parinfo[x0_ind+1].limits = [0.0, ylen+0.05]
;parinfo[x0_ind].fixed=1
;parinfo[x0_ind+1].fixed=1

if (keyword_set(freebulge)) then parinfo[bulge_ind].fixed = 0
if (keyword_set(freedisk)) then parinfo[2].fixed = 0

;default is to fix the sky
parinfo[skyInd].value = 0.0
if not keyword_set(free_sky) then $
  parinfo[skyInd].fixed = 1

start_params = parinfo[*].value

                                        ;normalize psf flux to one
data.psf /= total(data.psf)

                                        ;normalize maximum pixel value
                                        ;to maximum pixel flux
tot = max([total(data.image),10])
if( nparams lt 16) then $
  parinfo[0].value = tot/totalsersicflux(start_params) $
else begin
   rb=start_params[9]
   if rescale ne 0 then rb *= start_params[1]
    d = max([totalsersicflux([start_params[0:4],0.0,0.0,0.0]),1.0e-2])
    b = max([totalsersicflux([start_params[8], $
                              rb, $
                              start_params[10:12],0.0,0.0,0.0]),1.0e-2])

    parinfo[8].value = frac*tot / b
    parinfo[0].value = tot*(1.0-frac)/d
endelse
  start_params = parinfo[*].value
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
    ;print, start_params
    params=mpfit2dfun('pixelfluxpsf', x, y, data.image, $
                      data.ivar, $
                      start_params, parinfo=parinfo, $
                      functargs={psfImage:data.psf, $
                                 cutoff:1, rescale:rescale}, $
                      perror=errors, covar=covar, weights=data.ivar,$
                      dof=dof, bestnorm=bn, $;/quiet, $
                      xtol=1.0e-10, gtol=1.0e-8, ftol=1.0e-6, $
                      status=fitstat, maxiter=200, npegged=npg);, $
                      ;iterproc='iter_plot', ITERARGS={psf:data.psf}, yfit=fit)


    if( nparams ge 16 and times lt 4 ) then begin
        ;don't allow one set of parameters to be zero
        if (params[1] - parinfo[1].limits[0]) lt 1.e-8 then begin
            params[1] = start_params[1]*(1.-0.2*times)
            again = 1                
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
        if (rescale ne 0) and $
           (parinfo[9].limits[1] - params[9] lt 1.0e-5) then begin
           params[9] = (0.9-0.1*times)
           again=1
        endif
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
            params[8] = frac*tot/b
            again = 1
        endif
        ;if the fit is mostly bulge, but has a large, low level disk
        ; try making a small central exponential disk+larger dVc component
        ;if reset ne 1 and bulgetotot(params,cutoff=1) gt 0.6 and $
        ;  params[8] gt 1.05 and params[3] gt 0.3 and $
        ;  params[7] lt params[1] then begin
        ;    params = start_params
        ;    params[1] /= 8.0*_EXTRA.fracdev-3.0
        ;    d=max([totalsersicflux([1.0,params[1:4],0.0,0.0,0.0]),1.0e-2])
        ;    b=max([totalsersicflux([params[6:10],0.0,0.0,0.0]),1.0e-2])
        ;    params[0] = (1.0-frac)*tot/d
        ;    params[6] = frac*tot/b
        ;    again = 1
        ;    reset = 1
        ;endif
    endif 
    ;make q <=1 for disk and bulge
    if( (params[3] gt 1.0) and abs(params[0]) gt 1.0e-18 $
        and (params[1] gt 1.0e-12)) then begin
        params[3] = 1.0/params[3]
        params[phi_ind] += !pi*0.5
        params[1] /= params[3]
        if (nparams ge 16) and (rescale ne 0) then $
           params[9] *= params[3]
        if (again eq 1) then $
          params[3] = max([params[3], 1.1*parinfo[3].limits[0]])
        ;again  = 0
    endif 
    if( nparams ge 16 ) then begin
        if( params[11] gt 1.0  and abs(params[9]) gt 1.0e-18 $
            and params[8] gt 1.0e-12 ) then begin
            params[11] = 1.0/params[11]
            params[15] += !pi*0.5
            params[9] /= params[11]
            if (rescale ne 0) and (params[9] gt 1.005) then  begin
               again = 1
               params[9] = 0.8
            endif
            if (again eq 1) then $
              params[11] = max([params[11], 1.1*parinfo[11].limits[0]])
            ;again  = 0
        endif
        params[15] = params[15] mod !pi
    endif
    
    ;if the fit didn't work, try again
    if( fitstat eq 5 ) then begin
        again=1
        params[0] *= 0.9
        if( nparams ge 16) then begin
            params[8] *= 1.1
        endif
    endif
    params[phi_ind] = params[phi_ind] mod !pi
    temp = params
    times += 1
endwhile

if(times ge 5 and again eq 1) then begin
    print, 'reached max number of tries before fit was found'
    print, 'bulgeType: ', bulgeSersic
    print, 'status=', fitstat
    print, 'times=',times
    print, 'params=',params
    params[*] = -1.0
    fitstat = 0
endif


;put the rotation in units of pi
params[phi_ind] = params[phi_ind] mod !pi
if nparams ge 16 then params[15] = params[15] mod !pi


if nparams ge 16 then begin
    diskparams = [params[0:7]];,params[12:13],params[5]]
    bulgeparams = [params[8:15]];, params[12:13], params[11]]
endif else begin
    bulgeparams = [params[0:7]]
    diskparams = -1.0
 endelse
if keyword_set(rescale) then begin
   errors[9] = sqrt(diskparams[1]^2*errors[9]^2 + $
                    bulgeparams[1]^2*errors[1]^2 + $
                    2*diskparams[1]*bulgeparams[1]*covar[1,9])
   covar[*,9] *= diskparams[1] ;THIS IS WRONG!!!
   covar[9,*] *= diskparams[1] ;still WRONG
   bulgeparams[1] *= diskparams[1]
endif
;don't count null pixels in the final DOF to get the reduced chi squared
junk=where(data.ivar lt 1.0e-12, nullpix)
chisquare=bn/(dof - nullpix)
status=fitstat
degfree = dof-nullpix
skyVal = params[skyInd]

END
