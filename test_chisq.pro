;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; tests the chi_sq value of the high and low resolution 
; and high and low signal to noise images fits 
; against eachother
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO test_chisq
pref_set, 'idl_cpu_tpool_nthreads', 0, /commit
;the data that has no added noise
S0data = mrdfits('../output_01062011/keeplist39_sky/totalRAWSN_T00.00.fits',1)
;data with no added resolution
R0data = mrdfits('../output_01062011/keeplist39_sky/totalRAWRES00.00.fits',1)
;the original data
origdata = mrdfits('../output_01062011/keeplist39_sky/totalRAW_orig.fits',1)

;chi square of degraded model at high resolutions
chisq_ORIGres = dblarr(n_elements(S0data),5)
chisq_ORIGsn = dblarr(n_elements(R0data),5)

;chi square of high res model for degraded image
chisq_RESorig = dblarr(n_elements(S0data),5)
chisq_SNorig = dblarr(n_elements(R0data),5)

;chi squared of actual fits
chisq_RES = dblarr(n_elements(S0data), 5)
chisq_origr = dblarr(n_elements(S0data), 5)
chisq_SN = dblarr(n_elements(R0data), 5)
chisq_origsn = dblarr(n_elements(R0data), 5)

seed=systime(1)

for g=0, n_elements(origdata)-1 do begin

    res_match = where(S0data.GALID eq origdata[g].GALID)
    sn_match = where(R0data.GALID eq origdata[g].GALID)
    ;print, sn_match
    
    data = get_imivarpsf(origdata[g].RUN, origdata[g].CAMCOL, origdata[g].FIELD,$
                         origdata[g].ID, origdata[g].ROWC[2], origdata[g].COLC[2],$
                         rerun=origdata[g].RERUN)

    imsize = size(data.image, /dimensions)

    makegrid, imsize[0], imsize[1],x,y
    
    for gres=0,n_elements(res_match)-1 do begin
        gr = res_match[gres]
        add_fwhm = S0data[gr].ADD_FWHM
        chisq_RES[gr,*] = transpose([S0data[gr].CHISQ_BULGE, S0data[gr].CHISQ_DISK, $
                                     S0data[gr].CHISQ_DVC, S0data[gr].CHISQ_EXP, $
                                     S0data[gr].CHISQ_SERSIC])
        chisq_origr[gr,*]=transpose([origdata[g].CHISQ_BULGE, origdata[g].CHISQ_DISK, $
                                     origdata[g].CHISQ_DVC, origdata[g].CHISQ_EXP, $
                                     origdata[g].CHISQ_SERSIC])
        
        if add_fwhm lt 1.e-3 then continue
        ;print, add_fwhm, 'added fwhm'
        
        im_fft = fft(data.image, /double)
        noise = dblarr(imsize[0], imsize[1])
        noise[where(data.ivar ne 0)] = sqrt(1./data.ivar[where(data.ivar ne 0)])
        noise_fft = fft(noise, /double)
        psf_fft = fft(data.psf, /double)
        add_psf_correct = makepsf(51,51, (add_fwhm-1.)*origdata[g].PSF_FWHM[2])
                                        ;arrange the psf so the
                                        ;convolution comes out tin the
                                        ;correct quadrants
        add_psf=dblarr(imsize[0],imsize[1])
        add_psf[imsize[0]-51/2,imsize[1]-51/2] = $
          add_psf_correct[0:51/2-1,0:51/2-1]
        add_psf[0,imsize[1]-51/2] = add_psf_correct[51/2:51-1,$
                                                        0:51/2-1]
        add_psf[imsize[0]-51/2,0] = $
          add_psf_correct[0:51/2-1,51/2:51-1]
        add_psf[0,0] = add_psf_correct[51/2:51-1, 51/2:51-1]
        add_psf_fft = fft(add_psf, /double)

        new_image = (double(fft(im_fft*add_psf_fft, /double, /inverse)))*imsize[0]*imsize[1]
        new_psf = convolve(data.psf,add_psf_correct,/correlate) 
        noise = (double(fft(noise_fft*add_psf_fft, /double, $
                            /inverse)))*imsize[0]*imsize[1]
        new_ivar = dblarr(imsize[0],imsize[1])
        new_ivar[where(data.ivar ne 0)] = 1./(noise[where(data.ivar ne 0)]^2)
   
        chisq_ORIGres[gr,0] = total((data.image - pixelfluxpsf(x,y,S0data[gr].BULGEFIT,$
                                                               _EXTRA={psfImage:data.psf,$
                                                                       cutoff:1}))^2*data.ivar)/$
          origdata[g].DOF_BULGE
        chisq_ORIGres[gr,1] = total((data.image - pixelfluxpsf(x,y,S0data[gr].DISKFIT,$
                                                               _EXTRA={psfImage:data.psf,$
                                                                       cutoff:1}))^2*data.ivar)/$
          origdata[g].DOF_DISK
        chisq_ORIGres[gr,2] = total((data.image - pixelfluxpsf(x,y,S0data[gr].DVCFIT,$
                                                               _EXTRA={psfImage:data.psf,$
                                                                       cutoff:1}))^2*data.ivar)/$
          origdata[g].DOF_DVC
        chisq_ORIGres[gr,3] = total((data.image - pixelfluxpsf(x,y,S0data[gr].EXPFIT,$
                                                               _EXTRA={psfImage:data.psf,$
                                                                       cutoff:1}))^2*data.ivar)/$
          origdata[g].DOF_EXP
        chisq_ORIGres[gr,4] = total((data.image - pixelfluxpsf(x,y,S0data[gr].SERSICFIT,$
                                                               _EXTRA={psfImage:data.psf,$
                                                                       cutoff:1}))^2*data.ivar)/$
          origdata[g].DOF_SERSIC



        chisq_RESorig[gr,0] = total((new_image - pixelfluxpsf(x,y,origdata[g].BULGEFIT,$
                                                              _EXTRA={psfImage:new_psf,$
                                                                      cutoff:1}))^2*new_ivar)/$
          origdata[g].DOF_BULGE
        chisq_RESorig[gr,1] = total((new_image - pixelfluxpsf(x,y,origdata[g].DISKFIT,$
                                                              _EXTRA={psfImage:new_psf,$
                                                                      cutoff:1}))^2*new_ivar)/$
          origdata[g].DOF_DISK

        chisq_RESorig[gr,2] = total((new_image - pixelfluxpsf(x,y,origdata[g].DVCFIT,$
                                                              _EXTRA={psfImage:new_psf,$
                                                                      cutoff:1}))^2*new_ivar)/$
          origdata[g].DOF_DVC
        chisq_RESorig[gr,3] = total((new_image - pixelfluxpsf(x,y,origdata[g].EXPFIT,$
                                                              _EXTRA={psfImage:new_psf,$
                                                                      cutoff:1}))^2*new_ivar)/$
          origdata[g].DOF_EXP
        chisq_RESorig[gr,4] = total((new_image - pixelfluxpsf(x,y,origdata[g].SERSICFIT,$
                                                              _EXTRA={psfImage:new_psf,$
                                                                      cutoff:1}))^2*new_ivar)/$
          origdata[g].DOF_SERSIC
    endfor

    ;signal to noise case
    for gsn=0,n_elements(sn_match)-1 do begin
        gr = sn_match[gsn]
        add_noise = R0data[gr].ADD_NOISE
      
        chisq_SN[gr,*] = transpose([R0data[gr].CHISQ_BULGE[0], R0data[gr].CHISQ_DISK[0], $
                                  R0data[gr].CHISQ_DVC[0], R0data[gr].CHISQ_EXP[0], $
                                  R0data[gr].CHISQ_SERSIC[0]])
        chisq_origsn[gr,*]=transpose([origdata[g].CHISQ_BULGE, origdata[g].CHISQ_DISK, $
                                      origdata[g].CHISQ_DVC, origdata[g].CHISQ_EXP, $
                                      origdata[g].CHISQ_SERSIC])
        
        if add_noise lt 1.e-3 then continue
        ;print, add_noise, 'added noise'
        new_image = data.image
        new_image[where(data.ivar ne 0)] += ((randomn(seed, imsize[0], $
                                                      imsize[1]))[where(data.ivar ne 0)])*$
          sqrt(1./(data.ivar[where(data.ivar ne 0)]))*$
          sqrt(add_noise*add_noise-1.)
        new_ivar = data.ivar
        new_ivar[where(data.ivar ne 0)] = 1./((1./(data.ivar[where(data.ivar ne 0)]))*(add_noise*add_noise))

        chisq_ORIGsn[gr,0] = total((data.image - pixelfluxpsf(x,y,R0data[gr].BULGEFIT,$
                                                              _EXTRA={psfImage:data.psf,$
                                                                       cutoff:1}))^2*data.ivar)/$
          origdata[g].DOF_BULGE
        chisq_ORIGsn[gr,1] = total((data.image - pixelfluxpsf(x,y,R0data[gr].DISKFIT,$
                                                              _EXTRA={psfImage:data.psf,$
                                                                       cutoff:1}))^2*data.ivar)/$
          origdata[g].DOF_DISK
        chisq_ORIGsn[gr,2] = total((data.image - pixelfluxpsf(x,y,R0data[gr].DVCFIT,$
                                                              _EXTRA={psfImage:data.psf,$
                                                                      cutoff:1}))^2*data.ivar)/$
          origdata[g].DOF_DVC
        chisq_ORIGsn[gr,3] = total((data.image - pixelfluxpsf(x,y,R0data[gr].EXPFIT,$
                                                              _EXTRA={psfImage:data.psf,$
                                                                       cutoff:1}))^2*data.ivar)/$
          origdata[g].DOF_EXP
        chisq_ORIGsn[gr,4] = total((data.image - pixelfluxpsf(x,y,R0data[gr].SERSICFIT,$
                                                              _EXTRA={psfImage:data.psf,$
                                                                       cutoff:1}))^2*data.ivar)/$
          origdata[g].DOF_SERSIC
        
        chisq_SNorig[gr,0] = total((new_image - pixelfluxpsf(x,y,origdata[g].BULGEFIT,$
                                                             _EXTRA={psfImage:data.psf,$
                                                                     cutoff:1}))^2*new_ivar)/$
          origdata[g].DOF_BULGE
        chisq_SNorig[gr,1] = total((new_image - pixelfluxpsf(x,y,origdata[g].DISKFIT,$
                                                             _EXTRA={psfImage:data.psf,$
                                                                     cutoff:1}))^2*new_ivar)/$
          origdata[g].DOF_DISK
        chisq_SNorig[gr,2] = total((new_image - pixelfluxpsf(x,y,origdata[g].DVCFIT,$
                                                             _EXTRA={psfImage:data.psf,$
                                                                     cutoff:1}))^2*new_ivar)/$
          origdata[g].DOF_DVC
        chisq_SNorig[gr,3] = total((new_image - pixelfluxpsf(x,y,origdata[g].EXPFIT,$
                                                             _EXTRA={psfImage:data.psf,$
                                                                     cutoff:1}))^2*new_ivar)/$
          origdata[g].DOF_EXP
        chisq_SNorig[gr,4] = total((new_image - pixelfluxpsf(x,y,origdata[g].SERSICFIT,$
                                                             _EXTRA={psfImage:data.psf,$
                                                                     cutoff:1}))^2*new_ivar)/$
          origdata[g].DOF_SERSIC
        
    endfor

endfor


output_enRES = {GALID:'SDSSJhhmmss.ss+ddmmss.s', $
                ADD_FWHM:0.D, CHISQ_ADDFWHM:dblarr(5), $
                CHISQ_ORIG:dblarr(5), $
                CS_IMFWHM_ORIG:dblarr(5), $
                CS_IMORIG_FWHM:dblarr(5)}
output_RES = replicate( output_enRES, n_elements(S0data))
output_RES.GALID = S0data.GALID
output_RES.ADD_FWHM = S0data.ADD_FWHM
output_RES.CHISQ_ADDFWHM = transpose(chisq_RES)
output_RES.CHISQ_ORIG = transpose(chisq_origr)
output_RES.CS_IMFWHM_ORIG = transpose(chisq_RESorig)
output_RES.CS_IMORIG_FWHM = transpose(chisq_ORIGres)
mwrfits, output_RES, '../output_01062011/keeplist39_sky/RES_CScomp.fits', /create



output_enSN = {GALID:'SDSSJhhmmss.ss+ddmmss.s', $
                ADD_NOISE:0.D, CHISQ_ADDNOISE:dblarr(5), $
                CHISQ_ORIG:dblarr(5), $
                CS_IMNOISE_ORIG:dblarr(5), $
                CS_IMORIG_NOISE:dblarr(5)}
output_SN = replicate( output_enSN, n_elements(R0data))
output_SN.GALID = R0data.GALID
output_SN.ADD_NOISE = R0data.ADD_NOISE
output_SN.CHISQ_ADDNOISE = transpose(chisq_SN)
output_SN.CHISQ_ORIG = transpose(chisq_origsn)
output_SN.CS_IMNOISE_ORIG = transpose(chisq_SNorig)
output_SN.CS_IMORIG_NOISE = transpose(chisq_ORIGsn)
mwrfits, output_SN, '../output_01062011/keeplist39_sky/SN_CScomp.fits', /create


END
