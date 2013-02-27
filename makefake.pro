;____________________________________
;
; Claire Lackner
;  Jan 11 2010
;
;  makefake.pro
;
;  makes a fake image of a galaxy, either pure elliptical or
;  disk+bulge or disk+disk
;
; the arguments are disk, bulge, and params
;  if disk ne 0, then there should be disk
;  if bulge ne 0, there should a be sersic bulge
;  if bulge eq 0, there should be a disky bulge
;
; This gives the following for different galaxies
; Elliptical:  profile = 1
; Disk+Bulge:  profile = 2
; Disk+disk:   profile = 3
;
; the params is a structure consisting of disk:[8x1 array] and
;bulge:[8x1  array] of the input parameters, PSF gives the psf that
;the profiles should be convolved with, if PARAMS and/or PSF are not
;set, then randomly generated, reasonable parameters are used  
;
;____________________________________

FUNCTION posrand, mean, sig, seed
x = -1.0
while x LT 0 do begin
    x = randomn( seed, 1 )*sig + mean
endwhile
return, x
END

FUNCTION makefake, profile, xd, yd, in, PARAMS=PARAMS, PSF=PSF

;make the bulge
if( keyword_set(PARAMS) ) then begin
    if( profile eq 1 ) then inputs = PARAMS.bulge else $
      inputs=[[[PARAMS.disk]],[[PARAMS.bulge]]]
endif else begin
    seed = systime(1)
    ;input params need to be generated
    if( profile eq 1 or profile eq 2 ) then begin ;sersic bulge
        inputs = dblarr(8,1)
        inputs[0] = randomu(seed, 1)*10.0
        inputs[2] = (randomu(seed,1) + 0.6D)*3.2D
        inputs[1] = posrand( 12., 5., seed )
        inputs[3] = 1.0-randomu(seed,1)^2
        inputs[4] = 0.05*randomn(seed,1)
        inputs[5] = posrand(200.,3.,seed)
        inputs[6] = posrand(200.,3.,seed)
        inputs[7] = randomu(seed,1)*!PI*0.5
    endif
    if( profile eq 3 ) then begin ;disk pseudobulge
        inputs = dblarr(8,1)
        inputs[0] = randomu(seed)*10.
        inputs[1] = posrand(10.,5.,seed)
        inputs[2] = inputs[1]*posrand(0.10, 0.05, seed )
        inputs[3] = randomu(seed,1)*!PI*0.5
        inputs[4] = inputs[2] * 25.
        inputs[5:6] = abs(randomn(seed,2)*3.0+50.)
        inputs[7] = randomu(seed,1)*!PI*0.5
    endif
    xlen = fix(inputs[5]*2.)
    ylen = fix(inputs[6]*2.)
    if( profile eq 2 or profile eq 3) then begin ;outerdisk
        inputb = inputs
        inputs = dblarr( 8,2 )
        inputs[*,1] = inputb
        inputs[0,0] = inputb[0]*posrand(4.,0.6,seed)
        if( profile eq 2 ) then $
          inputs[1,0] = inputb[2]*(1.0+randomu(seed,1)*3.0) else $
          inputs[1,0] = inputb[1]*(1.0+randomu(seed,1)*3.0)
        inputs[2,0] = inputs[1,0]*posrand(0.10, 0.05, seed )
        if( profile eq 3 ) then $ ;keep inner disk inclined like outerdisk
          inputs[3,0] = inputs[3,1] + randomn(seed)*0.4 else $
          inputs[3,0] =  randomu(seed,1)*!PI*0.5
        inputs[4,0] = inputs[2,0] * 25.
        inputs[5:6,0] = inputb[5:6]
        if( profile eq 3 ) then $ ;keep inner disk inclined like outerdisk
          inputs[7,0] = inputs[7,1] + randomn(seed)*0.4 else $
          inputs[7,0] =  randomu(seed,1)*!PI*0.5
    endif
    
endelse


if(keyword_set(PSF) ) then fwhm=PSF else $
  fwhm = abs(randomn(seed)*0.5 + 2.6)


image = makegal( profile, inputs, psf=fwhm, imdimen=[xlen, ylen] )
xd = xlen
yd = ylen
in = inputs
return, image

END


