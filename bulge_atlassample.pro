;____________________________________
;
; Claire Lackner
;  Jan 14 2010
;
;  bulgeatlassample.pro
;
;  fits a bulge+disk to the atlas sample galaxies
;
; 
;____________________________________


;uses atlassample galaxies from atlasssampleDR7
PRO bulge_atlassample, filename, start, last, outputdir, imagedir, $
                       nosersic=nosersic, redo=redo, residuals=residuals
print, 'cpus: ',!CPU.TPOOL_NTHREADS     

gals = mrdfits(filename,1,range=[start,last-1])


;output entry for each object, the size of the array is for the bands,
;the default is five
if not keyword_set(redo) then begin
    output_entry = create_struct( gals[0], ['BULGEFIT', 'DISKFIT',$
                                            'SERSICFIT', 'CHISQ_BULGE', $
                                            'CHISQ_DISK', 'CHISQ_SERSIC', $
                                            'COVAR_BULGE', 'COVAR_DISK', $
                                            'COVAR_SERSIC', 'PERR_BULGE', $
                                            'PERR_DISK', 'PERR_SERSIC',$
                                            'MPFIT_STATUS', $
                                            'DOF_BULGE', 'DOF_DISK',$
                                            'DOF_SERSIC', 'DOF_DVC',$
                                            'DOF_EXP',$
                                            'EXPFIT', 'DVCFIT', $
                                            'CHISQ_EXP', 'CHISQ_DVC', $
                                            'PERR_EXP', 'PERR_DVC', $
                                            'COVAR_EXP', 'COVAR_DVC',$
                                            'FRACDEV',$
                                            'XCROP', 'YCROP',$
                                            'XLEN', 'YLEN',$
                                            'DVC_BTT', 'EXP_BTT',$
                                            'MAD_SKY', 'MAD_SERSIC', $
                                            'MAD_SERSIC_MASK', $
                                            'MAD_DVCB', $
                                            'MAD_DVCB_MASK',$
                                            'MAD_EXPB',$
                                            'MAD_EXPB_MASK',$
                                            'MAD_EXP',$
                                            'MAD_EXP_MASK',$
                                            'MAD_DVC',$
                                            'MAD_DVC_MASK',$
                                            'CHISQ_BULGE_MASK',$
                                            'CHISQ_DISK_MASK',$
                                            'CHISQ_EXP_MASK',$
                                            'CHISQ_SERSIC_MASK',$
                                            'CHISQ_DVC_MASK',$
                                            'DOF_BULGE_MASK',$
                                            'DOF_DISK_MASK',$
                                            'DOF_EXP_MASK',$
                                            'DOF_SERSIC_MASK',$
                                            'DOF_DVC_MASK', $
                                            'SN_REFF_SERSIC', $
                                            'SKY_BULGE', $
                                            'SKY_BULGE_ERR', $
                                            'SKY_BULGE_COVAR', $
                                            'DVC_BTT_ERR', $
                                            'EXP_BTT_ERR' ], $
                                  dblarr(16), dblarr(16), $
                                  dblarr(8), 0.0D, $
                                  0.0D, 0.0D, dblarr(16,16), $
                                  dblarr(16,16), dblarr(8,8), $
                                  dblarr(16), dblarr(16), $
                                  dblarr(8), dblarr(5),$
                                  0.0D,0.0D,$
                                  0.0D,0.0D,0.0D,$
                                  dblarr(8), dblarr(8),$
                                  0.0D, 0.0D, $
                                  dblarr(8),dblarr(8),$
                                  dblarr(8,8), dblarr(8,8), 0.0D, $
                                  0L, 0L, 0L, 0L, 0.0D, 0.0D,$
                                  0.0D, 0.0D, 0.0D, 0.0D, 0.0D, $
                                  0.0D, 0.0D, 0.0D, 0.0D, 0.0D, $
                                  0.0D, 0.0D, 0.0D, 0.0D, 0.0D, 0.0D,$
                                  0.0D, 0.0D, 0.0D, 0.0D, 0.0D, 0.0D,$
                                  0.0D, 0.0D, dblarr(17), 0.0D, 0.0D )
endif else output_entry = create_struct(gals[0])

;insert known values into output strucutre array
output = replicate( output_entry, n_elements(gals) )
struct_assign, gals, output
;output.BANDS = (where(notignore ne 0))#(intarr(n_elements(gals))+1)

;file for writing out times
;openw, 1, 'timing.dat'

;loop over each image
for i=0L, n_elements(gals)-1L do begin

    print, gals[i].IDENT, ': ', i+1, ' of ', n_elements(gals)
    t1=systime(1)
        
    data = get_imivarpsf(gals[i].IDENT, strtrim(gals[i].FILENAME), imagedir)

    origsize = size(data.image, /dimensions )
    if n_elements(origsize) eq 1 then continue
    xcrop = 0
    ycrop = 0
    origdata=data
    data = cropimage_2(data.image, data.ivar, data.psf, $
                     data.mask, $
                     x0=xcrop, y0=ycrop)
    ;x0 = 0
    ;y0 = 0
    imsize=size(data.image, /dimensions)
    if n_elements(imsize) lt 2 then begin
       print, '   cropped size: 0, skipping'
       continue
    endif

    print, '   cropped size:', imsize[0],'x',imsize[1]
    output[i].XLEN=imsize[0]
    output[i].YLEN=imsize[1]
    output[i].XCROP=xcrop
    output[i].YCROP=ycrop
        
    makegrid, imsize[0], imsize[1], x, y

    if( not keyword_set(nosersic) ) then begin
        diskbulgefit, nod, params, data.image, data.psf, data.ivar, $
          cs, covar, errors, fitstat, $
          dof, sky, bulgeSersic=1.0, $ ;4.0*fracDev + 1.0, $ 
          nodisk=1, freebulge=1;, /free_sky
                   
        params[5] += xcrop
        params[6] += ycrop
        output[i].SERSICFIT[*] = params
        output[i].CHISQ_SERSIC = cs
        output[i].PERR_SERSIC[0:7] = errors[0:7]
        output[i].COVAR_SERSIC[0:7,0:7] = covar[0:7,0:7]          
        output[i].MPFIT_STATUS[2] = fitstat
        output[i].DOF_SERSIC = dof
        ;output[i].SKY_SERSIC = sky
        ;output[i].SKY_SERSIC_ERR = errors[8]
        ;output[i].SKY_SERSIC_COVAR = covar[8,0:8]
    endif

;;                                         ;fit dVc profile
;;     diskbulgefit,nod,params,data.image,data.psf,data.ivar,$
;;       cs,covar,errors,fitstat,dof, sky, bulgeSersic=4.0, $
;;       nodisk=1

;;     dvcparams = params
;;     params[5] += xcrop
;;     params[6] += ycrop
;;     output[i].DVCFIT = params
;;     output[i].CHISQ_DVC = cs
;;     output[i].PERR_DVC[0:7] = errors[0:7]
;;     output[i].COVAR_DVC[0:7,0:7] = covar[0:7,0:7]          
;;     output[i].MPFIT_STATUS[3] = fitstat
;;     output[i].DOF_DVC = dof

        
;;     diskbulgefit,nod,params,data.image,data.psf,data.ivar,$
;;       cs,covar,errors,fitstat,dof, sky, bulgeSersic=1.0, $
;;       nodisk=1

;;     expparams = params
;;     params[5] += xcrop
;;     params[6] += ycrop
;;     output[i].EXPFIT[*] = params
;;     output[i].CHISQ_EXP = cs
;;     output[i].PERR_EXP[0:7] = errors[0:7]
;;     output[i].COVAR_EXP[0:7,0:7] = covar[0:7,0:7]          
;;     output[i].MPFIT_STATUS[4] = fitstat
;;     output[i].DOF_EXP = dof


;;    fracDev = getfracdev(data.image, data.psf,data.ivar, $
;;                          expparams, dvcparams)
;;    output[i].fracdev = fracdev
;; ;if fracdev doesn't work, replace it with 0.5
;;    if fracdev ne fracdev then fracdev = 0.5


;; ;    guess = [output[i].EXPFIT[0:4], output[i].EXPFIT[7], $
;; ;             output[i].DVCFIT[0:4], output[i].DVCFIT[7], $
;; ;             output[i].DVCFIT[5:6]-[xcrop,ycrop]]
;; ;    guess[6] *= (0.1+fracDev*0.8)
;; ;    guess[0] *= (1.-(0.1+fracDev*0.8))
 
;;     diskbulgefit, dbparams, bbparams, $
;;       data.image, data.psf,$
;;       data.ivar, chsqdb, covardb, errdb, stat, dofdb, sky, bulgeSersic=4.0D, $
;;       /rescale, $
;;       ;/free_coords, $ ;different centroids
;;       ;/free_sky, $ ;fit for the sky
;;       _EXTRA={Reff:(1.-fracDev)*(output[i].EXPFIT[1]) + $
;;                     fracDev*(output[i].DVCFIT[1]), $
;;               fracdev:fracDev, $
;;               q:output[i].DVCFIT[3], $
;;               phi:output[i].EXPFIT[7]}
;;     output[i].MPFIT_STATUS[0] = stat
;;     output[i].SKY_BULGE = sky
;;     output[i].SKY_BULGE_ERR = errdb[16]
;;     output[i].SKY_BULGE_COVAR = covardb[16,0:16]

;;  ;       print, 'timed: ',systime(1)-t1
;; ;    guess[8]=1.0
;;     diskbulgefit, doparams, diparams, $
;;       data.image, data.psf, $
;;       data.ivar, chsqdd, covardd, errdd, stat, dofdd, sky, bulgeSersic=1.0D, $
;;       /rescale, $
;;       ;/free_coords, $                   ;different centroids
;;       _EXTRA={Reff:output[i].EXPFIT[1], $
;;               q:output[i].EXPFIT[3], $
;;               phi:output[i].EXPFIT[7]}
;;     output[i].MPFIT_STATUS[1] = stat
    
;;     ;correct for rescaling the bulge reff to the disk reff
;;     ;THIS is a HACK since I haven't dealt with the errors
;;     ;or the covariances!!!!
;;     ;diparams[1] *= doparams[1]
;;     ;bbparams[1] *= dbparams[1]

;;     doparams[5:6] += [xcrop, ycrop]
;;     diparams[5:6] += [xcrop, ycrop]
;;     if doparams[1] lt diparams[1] then begin
;;        tempdd = doparams
;;        doparams = diparams
;;        diparams = tempdd
       
;;        tempdd = errdd[0:7]
;;        errdd[0:7] = errdd[8:15]
;;        errdd[8:15 ] = tempdd

;;        tempdd = covardd
;;        tempdd[0:7,*] = covardd[8:15,*]
;;        tempdd[8:15,*] = covardd[0:7,*]
;;        tempdd2 = tempdd
;;        tempdd[*,0:7] = tempdd2[*,8:15]
;;        tempdd[*,8:15] = tempdd2[*,0:7]
;;        covardd = tempdd
;;     endif

;;     output[i].DISKFIT[*] = [doparams[*],$
;;                             diparams[*] ];[doparams[0:4],doparams[7],$
;;                                         ;   diparams[0:4],diparams[7],$
;;                                         ;   diparams[5:6]+[xcrop,ycrop]]
;;     dbparams[5:6] += [xcrop, ycrop]
;;     bbparams[5:6] += [xcrop, ycrop]
;;     output[i].BULGEFIT[*] = [dbparams[*], $
;;                              bbparams[*]] ;[dbparams[0:4],dbparams[7],$
;;                                         ;   bbparams[0:4],bbparams[7],$
;;                                         ;  bbparams[5:6]+[xcrop,ycrop]]
;;     output[i].PERR_BULGE[0:15] = errdb[0:15]
;;     output[i].PERR_DISK[0:15] = errdd[0:15]
;;     output[i].COVAR_DISK[0:15,0:15] = covardd[0:15,0:15]
;;     output[i].COVAR_BULGE[0:15,0:15] = covardb[0:15,0:15]
;;     output[i].CHISQ_BULGE = chsqdb
;;     output[i].CHISQ_DISK = chsqdd
;;     output[i].DOF_BULGE = dofdb
;;     output[i].DOF_DISK = dofdd
;;     output[i].DVC_BTT = bulgetotot(output[i].BULGEFIT, cutoff=1)
;;     output[i].EXP_BTT = bulgetotot(output[i].DISKFIT, cutoff=1)
;;     output[i].DVC_BTT_ERR = bulgetotot_error(output[i].BULGEFIT,$
;;                                              output[i].COVAR_BULGE)
;;     output[i].EXP_BTT_ERR = bulgetotot_error(output[i].DISKFIT,$
;;                                              output[i].COVAR_DISK)

    
    
;;     keep = where(origdata.mask,complement=notkeep)
;;     output[i].MAD_SKY = median(abs(origdata.image[notkeep]))
;;     output[i].DOF_BULGE_MASK = n_elements(keep) - $
;;       (n_elements(data.image) - output[i].DOF_BULGE)
;;     output[i].DOF_DISK_MASK = n_elements(keep) - $
;;       (n_elements(data.image) - output[i].DOF_DISK)
;;     output[i].DOF_SERSIC_MASK = n_elements(keep) - $
;;       (n_elements(data.image) - output[i].DOF_SERSIC)
;;     output[i].DOF_DVC_MASK = n_elements(keep) - $
;;       (n_elements(data.image) - output[i].DOF_DVC)
;;     output[i].DOF_EXP_MASK = n_elements(keep) - $
;;       (n_elements(data.image) - output[i].DOF_EXP)

;;     ;put the model images in a file
;;     if keyword_set(residuals) then begin
;;         modName = string( output[i].GALID, format='("models/M",i09,".fits")')
;;         modName = outputdir+modName
;;     endif
;;     makegrid, origsize[0], origsize[1], x, y
        
;;     model = pixelfluxpsf(x,y,$
;;                          [output[i].BULGEFIT,output[i].SKY_BULGE], $
;;                          _EXTRA={cutoff:1,psfImage:data.psf})
    
;;     output[i].MAD_DVCB = median(abs(origdata.image-model))
;;     output[i].MAD_DVCB_MASK = median(abs((origdata.image-model)[keep]))
;;     output[i].CHISQ_BULGE_MASK = total((origdata.ivar[keep])* $
;;                                        (origdata.image[keep]-model[keep])^2)/$
;;       output[i].DOF_BULGE_MASK

;;      ;put the model images in a file
;;     if keyword_set(residuals) then begin
;;         mwrfits, model, modName, /create
;;     endif
;;     model = pixelfluxpsf(x,y,$
;;                          [output[i].DISKFIT,0.0], $
;;                          _EXTRA={cutoff:1,psfImage:data.psf})
;;     output[i].MAD_EXPB = median(abs(origdata.image-model))
;;     output[i].MAD_EXPB_MASK = median(abs((origdata.image-model)[keep]))
;;     output[i].CHISQ_DISK_MASK = total((origdata.ivar[keep])* $
;;                                       (origdata.image[keep]-model[keep])^2)/$
;;       output[i].DOF_DISK_MASK

;;      ;put the model images in a file
;;     if keyword_set(residuals) then begin
;;         mwrfits, model, modName
;;     endif
;;     model = pixelfluxpsf(x,y,$
;;                          [output[i].DVCFIT,0.0], $
;;                          _EXTRA={cutoff:1,psfImage:data.psf})
;;     output[i].MAD_DVC = median(abs(origdata.image-model))
;;     output[i].MAD_DVC_MASK = median(abs((origdata.image-model)[keep]))
;;     output[i].CHISQ_DVC_MASK = total((origdata.ivar[keep])* $
;;                                      (origdata.image[keep]-model[keep])^2)/$
;;       output[i].DOF_DVC_MASK

;;      ;put the model images in a file
;;     if keyword_set(residuals) then begin
;;         mwrfits, model, modName
;;     endif
;;     model = pixelfluxpsf(x,y,$
;;                          [output[i].EXPFIT,0.0], $
;;                          _EXTRA={cutoff:1,psfImage:data.psf})
;;     output[i].MAD_EXP = median(abs(origdata.image-model))
;;     output[i].MAD_EXP_MASK = median(abs((origdata.image-model)[keep]))
;;     output[i].CHISQ_EXP_MASK = total((origdata.ivar[keep])* $
;;                                      (origdata.image[keep]-model[keep])^2)/$
;;       output[i].DOF_EXP_MASK

;;      ;put the model images in a file
;;     if keyword_set(residuals) then begin
;;         mwrfits, model, modName
;;     endif
;;     model = pixelfluxpsf(x,y,$
;;                          [output[i].SERSICFIT,0.0], $
;;                          _EXTRA={cutoff:1,psfImage:data.psf})
;;     output[i].MAD_SERSIC = median(abs(origdata.image-model))
;;     output[i].MAD_SERSIC_MASK = median(abs((origdata.image-model)[keep]))
;;     output[i].CHISQ_SERSIC_MASK = total((origdata.ivar[keep])* $
;;                                         (origdata.image[keep]-model[keep])^2)/$
;;       output[i].DOF_SERSIC_MASK
    
;;                                         ;get the S/N with Reff for the Sersic profile
;;     dist_ellipse, ellip, [origsize[0], origsize[1]], $
;;       output[i].SERSICFIT[5], output[i].SERSICFIT[6], $
;;       1./output[i].SERSICFIT[3], $
;;       output[i].SERSICFIT[7]*180.0/!pi-90.0

;;     inside=where(ellip le $
;;                  output[i].SERSICFIT[1], ins)
;;     if(ins gt 0) then begin
;;         signal = total(origdata.image[inside])
;;         noise = sqrt(total(1.0/origdata.ivar[inside]))
;;         output[i].SN_REFF_SERSIC = signal/noise
;;     endif
;;      ;put the model images in a file
;;     if keyword_set(residuals) then begin
;;         mwrfits, model, modName
;;     endif

endfor
;close,1
;write output fits file: name of file=
append = string( start, last-1, format='("RAWFIT",i05,".",i05,".fits")')
outname = outputdir+append


mwrfits, output, outname, /create
print, 'writing output file: '+outname

END
