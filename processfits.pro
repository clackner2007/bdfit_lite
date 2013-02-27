;________________________________________
; Claire Lackner
; Mar 17 2010
;
; processes output fits from bulge_atlassample.pro
;
; creats new tables for fit parameters, as well as images+residuals
; for the fits
;
;---------------------------------------------


FUNCTION getChiSquare, image, ivar, model_params, _EXTRA=_EXTRA

makegrid, (size(image, /dimensions))[0], (size(image, /dimensions))[1], x, y
modelimage = pixelfluxpsf(x,y,model_params,_EXTRA=_EXTRA)
return, total((image-modelimage)^2*ivar)

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO processfits, inputfolder, inputfiles, outputRAW, outputtable, images=images, $
                 allbands=allbands, no_colors=no_colors

;outputRAW = inputfolder+'/totalRAW.fits'
;pref_set, 'idl_cpu_tpool_nthreads', 0, /commit

if( not file_test(outputRAW) ) then concatdata, inputfolder, inputfiles, outputRAW

print, outputRAW
data = mrdfits( outputRAW, 1 )

nfits=5
if(keyword_set(allbands)) then begin
    nbands = 5 ;number of bands in output, all of them
    fitbands = data[0].BANDS      ;where those bands are in the output
    rband = where(fitbands eq 2 )
;    fitbands = [1,2,3]
endif else begin
    nbands = n_elements(data[0].BANDS)
    rband = where(data[0].BANDS eq 2) ;change for color runs
    fitbands = indgen(nbands)
endelse
    
new_entry = {GALID:data[0].GALID, Z:0.0D, $
             RA:0.0D, DEC:0.0D, RUN:1L, CAMCOL:1L, FIELD:1L, $
             ID:1L, RERUN:data[0].RERUN, BULGE_SERS_N:dblarr(nfits,nbands), $
             BULGE_REFF:dblarr(nfits,nbands), BULGE_REFF_ERR:dblarr(nfits,nbands), $
             BULGE_TO_TOT:dblarr(nfits,nbands), BULGE_TO_TOT_ERR:dblarr(nfits,nbands), $
             BULGE_Q:dblarr(nfits,nbands), BULGE_Q_ERR:dblarr(nfits,nbands), $
             DISK_SCALE:dblarr(nfits,nbands), DISK_SCALE_ERR:dblarr(nfits,nbands), $
             DISK_INC:dblarr(nfits,nbands), DISK_INC_ERR:dblarr(nfits,nbands),$
             DISK_IO:dblarr(nfits,nbands), BULGE_IO:dblarr(nfits,nbands), $
             DISK_IO_ERR:dblarr(nfits,nbands), $
             SERSIC_N_ERR:dblarr(nbands),$
             BULGE_IO_ERR:dblarr(nfits,nbands), $
             BULGE_FLUX_ERR:dblarr(nfits,5), $ 
             DISK_FLUX_ERR:dblarr(nfits,5), $ ;in counts
             TOTAL_FLUX_ERR:dblarr(nfits,5), $ ;in counts
             BULGE_FLUX_COUNTS:dblarr(nfits,5), $ ;in counts
             DISK_FLUX_COUNTS:dblarr(nfits,5), $  ;in counts
             TOTAL_FLUX_COUNTS:dblarr(nfits,5), $  ;in counts
             DISK_PHI:dblarr(nfits,nbands), BULGE_PHI:dblarr(nfits,nbands),$
             CEN_X:dblarr(nfits,nbands), CEN_Y:dblarr(nfits,nbands), $
             CHI_SQ:dblarr(nfits,nbands), BULGE_REL_RBAND:dblarr(nfits,5), $
             DISK_REL_RBAND:dblarr(nfits,5), $
             BULGE_MAG:dblarr(nfits,5), DISK_MAG:dblarr(nfits,5), $
             TOTAL_MAG:dblarr(nfits,5), $
             EXTINCTION:dblarr(5), $ ;copy of SDSS DR7 galactic extinction
             PIXSCALE:dblarr(5), NMGYPERCOUNT:dblarr(5), $
             PHI_OFFSET_RAD:dblarr(5), ROWC:dblarr(5),$
             COLC:dblarr(5), PSF_FWHM:dblarr(5),$
             MPFIT_STATUS:intarr(nfits,nbands), BANDS:data[0].BANDS, $
             DOF:dblarr(nfits,nbands), SDSS_CMODEL_MAG:dblarr(5),$
             SDSS_AB_EXP:dblarr(5), SDSS_AB_DEV:dblarr(5), $
             SDSS_PETRO_FLUX:dblarr(5), SDSS_PETRO_FLUX_ERR:dblarr(5),$
             SKY_COUNTS:dblarr(5), SKY_ERR:dblarr(5), GAIN:dblarr(5), $
             ATLAS_X:0,ATLAS_Y:0, CHI_SQ_SCALED:dblarr(nfits,5), $
             DOF_SCALED:dblarr(nfits,5),BULGE_CENPIX:dblarr(nfits,5),$
             DISK_CENPIX:dblarr(nfits,5),B_MINUS_D_MAX:dblarr(nfits,5),$
             B_MINUS_D_MIN:dblarr(nfits,5)}


table = replicate(new_entry, n_elements(data) )
for ii=0L, n_elements(data)-1 do begin
    table[ii].GALID = sdss_galid(data[ii].RA, data[ii].DEC)
endfor
;table.GALID = data.GALID
table.Z = data.Z
table.RA = data.RA
table.DEC = data.DEC
table.RUN = data.RUN
table.CAMCOL = data.CAMCOL
table.FIELD = data.FIELD
table.ID = data.ID
table.RERUN = data.RERUN

fb=n_elements(fitbands)                 ;how many bands were fit, not
                                        ;the same as how many in the output
if( fb ne 1) then begin
    table.CHI_SQ[*,fitbands,*] = [reform(data.CHISQ_BULGE,1,fb,n_elements(data)),$
                                  reform(data.CHISQ_DISK,1,fb,n_elements(data)),$
                                  reform(data.CHISQ_DVC,1,fb,n_elements(data)),$
                                  reform(data.CHISQ_EXP,1,fb,n_elements(data)),$
                                  reform(data.CHISQ_SERSIC,1,fb,n_elements(data))]
    if( tag_exist(data, "DOF_BULGE") ) then $ 
      table.DOF[*,fitbands,*] = [reform(data.DOF_BULGE,1,fb,n_elements(data)),$
                                 reform(data.DOF_DISK,1,fb,n_elements(data)),$
                                 reform(data.DOF_DVC,1,fb,n_elements(data)),$
                                 reform(data.DOF_EXP,1,fb,n_elements(data)),$
                                 reform(data.DOF_SERSIC,1,fb,n_elements(data))]
endif else begin
    table.CHI_SQ[*,fitbands,*] = [reform(data.CHISQ_BULGE,1,n_elements(data)),$
                                  reform(data.CHISQ_DISK,1,n_elements(data)),$
                                  reform(data.CHISQ_DVC,1,n_elements(data)),$
                                  reform(data.CHISQ_EXP,1,n_elements(data)),$
                                  reform(data.CHISQ_SERSIC,1,n_elements(data))]
    if( tag_exist(data, "DOF_BULGE") ) then $ 
      table.DOF[*,fitbands,*] = [reform(data.DOF_BULGE,1,n_elements(data)),$
                                 reform(data.DOF_DISK,1,n_elements(data)),$
                                 reform(data.DOF_DVC,1,n_elements(data)),$
                                 reform(data.DOF_EXP,1,n_elements(data)),$
                                 reform(data.DOF_SERSIC,1,n_elements(data))]
endelse
table.MPFIT_STATUS[*,fitbands,*] = data.MPFIT_STATUS

;n=4 fit
table.BULGE_SERS_N[0,fitbands] = data.BULGEFIT[8,*]
table.BULGE_REFF[0,fitbands] = data.BULGEFIT[7,*]
table.BULGE_REFF_ERR[0,fitbands] = data.PERR_BULGE[7,*]
table.BULGE_Q[0,fitbands] = data.BULGEFIT[9,*]
table.BULGE_Q_ERR[0,fitbands] = data.PERR_BULGE[9,*]
table.DISK_SCALE[0,fitbands] = data.BULGEFIT[1,*]*0.59684 ; exponential scale length
table.DISK_SCALE_ERR[0,fitbands] = data.PERR_BULGE[1,*]*0.59684 ;exp scale length
table.DISK_INC[0,fitbands] = acos(data.BULGEFIT[3,*])
table.DISK_INC_ERR[0,fitbands] = data.PERR_BULGE[3,*]/(sqrt(1.-data.BULGEFIT[3,*]^2))
table.DISK_IO[0,fitbands] = data.BULGEFIT[0,*]*exp(kappa_sersic(1.0))
table.BULGE_IO[0,fitbands] =data.BULGEFIT[6,*]*exp(kappa_sersic(data.BULGEFIT[8,*]))
table.DISK_IO_ERR[0,fitbands] =data.PERR_BULGE[0,*]*exp(kappa_sersic(1.0))
table.BULGE_IO_ERR[0,fitbands] =(data.PERR_BULGE[6,*]*$
                          exp(kappa_sersic(data.BULGEFIT[8,*])))
table.BULGE_PHI[0,fitbands] = data.BULGEFIT[11,*] 
table.DISK_PHI[0,fitbands] = data.BULGEFIT[5,*]
table.CEN_X[0,fitbands] = data.BULGEFIT[12,*]
table.CEN_Y[0,fitbands] = data.BULGEFIT[13,*]

;no error on the disk inclinations
;n=1 fit
table.BULGE_SERS_N[1,fitbands] = data.DISKFIT[8,*]
table.BULGE_REFF[1,fitbands] = data.DISKFIT[7,*]
table.BULGE_REFF_ERR[1,fitbands] = data.PERR_DISK[7,*]
table.BULGE_Q[1,fitbands] = data.DISKFIT[9,*]
table.BULGE_Q_ERR[1,fitbands] = data.PERR_DISK[9,*]
table.DISK_SCALE[1,fitbands] = data.DISKFIT[1,*]*0.59684 ; exponential scale length
table.DISK_SCALE_ERR[1,fitbands] = data.PERR_DISK[1,*]*0.59684 ;exp scale length
table.DISK_INC[1,fitbands] = acos(data.DISKFIT[3,*])
table.DISK_INC_ERR[1,fitbands] = data.PERR_DISK[3,*]/(sqrt(1.-data.DISKFIT[3,*]^2))
table.DISK_IO[1,fitbands] = data.DISKFIT[0,*]*exp(kappa_sersic(1.0))
table.BULGE_IO[1,fitbands] =data.DISKFIT[6,*]*exp(kappa_sersic(data.DISKFIT[8,*]))
table.DISK_IO_ERR[1,fitbands] =data.PERR_DISK[0,*]*exp(kappa_sersic(1.0))
table.BULGE_IO_ERR[1,fitbands] =(data.PERR_DISK[6,*]*$
                          exp(kappa_sersic(data.DISKFIT[8,*])))
table.BULGE_PHI[1,fitbands] = data.DISKFIT[11,*] 
table.DISK_PHI[1,fitbands] = data.DISKFIT[5,*]
table.CEN_X[1,fitbands] = data.DISKFIT[12,*]
table.CEN_Y[1,fitbands] = data.DISKFIT[13,*]

;dVc only
table.BULGE_SERS_N[2,fitbands] = data.DVCFIT[2,*]
table.BULGE_REFF[2,fitbands] = data.DVCFIT[1,*]
table.BULGE_REFF_ERR[2,fitbands] = data.PERR_DVC[1,*]
table.BULGE_Q[2,fitbands] = data.DVCFIT[3,*]
table.BULGE_Q_ERR[2,fitbands] = data.PERR_DVC[3,*]
table.BULGE_IO[2,fitbands] =data.DVCFIT[0,*]*exp(kappa_sersic(4.0))
table.BULGE_IO_ERR[2,fitbands] =(data.PERR_DVC[0,*]*$
                          exp(kappa_sersic(4.0)))
table.BULGE_PHI[2,fitbands] = data.DVCFIT[7,*] 
table.CEN_X[2,fitbands] = data.DVCFIT[5,*]
table.CEN_Y[2,fitbands] = data.DVCFIT[6,*]

;exp disk only
table.DISK_SCALE[3,fitbands] = data.EXPFIT[1,*]*0.59684 ; exponential scale length
table.DISK_SCALE_ERR[3,fitbands] = data.PERR_EXP[1,*]*0.59684 ;exp scale length
table.DISK_INC[3,fitbands] = acos(data.EXPFIT[3,*])
table.DISK_INC_ERR[3,fitbands] = data.PERR_EXP[3,*]/(sqrt(1.-data.EXPFIT[3,*]^2))
table.DISK_IO[3,fitbands] = data.EXPFIT[0,*]*exp(kappa_sersic(1.0))
table.DISK_IO_ERR[3,fitbands] =data.PERR_EXP[0,*]*exp(kappa_sersic(1.0))
table.DISK_PHI[3,fitbands] = data.EXPFIT[7,*]
table.CEN_X[3,fitbands] = data.EXPFIT[5,*]
table.CEN_Y[3,fitbands] = data.EXPFIT[6,*]


;sersic fit only
table.SERSIC_N_ERR[fitbands] = data.PERR_SERSIC[2,*]
table.BULGE_SERS_N[4,fitbands] = data.SERSICFIT[2,*]
table.BULGE_REFF[4,fitbands] = data.SERSICFIT[1,*]
table.BULGE_REFF_ERR[4,fitbands] = data.PERR_SERSIC[1,*]
table.BULGE_Q[4,fitbands] = data.SERSICFIT[3,*]
table.BULGE_Q_ERR[4,fitbands] = data.PERR_SERSIC[3,*]
table.BULGE_IO[4,fitbands] =data.SERSICFIT[0,*]*exp(kappa_sersic(data.SERSICFIT[2,*]))
table.BULGE_IO_ERR[4,fitbands] =(data.PERR_SERSIC[0,*]*$
                                 exp(kappa_sersic(data.SERSICFIT[2,*])))
table.BULGE_PHI[4,fitbands] = data.SERSICFIT[7,*] 
table.CEN_X[4,fitbands] = data.SERSICFIT[5,*]
table.CEN_Y[4,fitbands] = data.SERSICFIT[6,*]


;get the extinctions and the nanomaggies per count from the objects
; extras = sdss_readobjlist(data.RUN, data.CAMCOL, $
;                           data.FIELD, data.ID, $
;                           rerun=data.RERUN, /silent, $
;                           _EXTRA={select_tags:['EXTINCTION',$
;                                                'NMGYPERCOUNT', 'PHI_OFFSET']})
;                                         ;carry the extinction on
; table.EXTINCTION = extras.EXTINCTION
; table.NMGYPERCOUNT = extras.NMGYPERCOUNT
; table.PHI_OFFSET_RAD = extras.PHI_OFFSET*!pi/180.0
table.EXTINCTION = data.EXTINCTION
table.NMGYPERCOUNT =data.NMGYPERCOUNT
table.PHI_OFFSET_RAD = data.PHI_OFFSET_RAD
table.PIXSCALE = data.PIXSCALE
table.ROWC = data.ROWC
table.COLC = data.COLC
table.PSF_FWHM = data.PSF_FWHM
table.SDSS_CMODEL_MAG = sdss_asinh_mag(data.CMODELFLUX*1.0e-9,filter=[0,1,2,3,4])
table.SDSS_AB_EXP = data.AB_EXP
table.SDSS_AB_DEV = data.AB_DEV
table.SDSS_PETRO_FLUX = data.PETROFLUX
table.SDSS_PETRO_FLUX_ERR = 1./sqrt(data.PETROFLUX_IVAR)



bands=['u','g','r','i','z']
fiducial = 'r'
;rband = where(fitbands eq 2)

for i=0L, n_elements(data)-1 do begin

    ;get the flux error(assume only rband)
    ;take a scaled r-band model and do a linear fit for each of the
    ;other bands, and the ncompute the magnitudes in these bands,
    ;currently, these are Pogson magnitudes, these have NOT been
    ;extinction corrected for the galaxy
    if not keyword_set(no_colors) then begin
        for color=0L, n_elements(bands)-1 do begin
                                        ;don't do the rband
            if bands[color] eq fiducial then begin
                table[i].BULGE_REL_RBAND[*,color] = 1.0
                table[i].DISK_REL_RBAND[*,color] = 1.0
                
                
                table[i].BULGE_FLUX_ERR[0,color]  = totalsersicflux_error(data[i].BULGEFIT,$
                                                                          data[i].COVAR_BULGE,$
                                                                          /bulge, cutoff=1)
                table[i].BULGE_FLUX_ERR[1,color] = totalsersicflux_error(data[i].DISKFIT,$
                                                                         data[i].COVAR_DISK,$
                                                                         /bulge, cutoff=1)
                table[i].DISK_FLUX_ERR[0,color]  = totalsersicflux_error(data[i].BULGEFIT,$
                                                                         data[i].COVAR_BULGE,$
                                                                         /disk, cutoff=1)
                table[i].DISK_FLUX_ERR[1,color] = totalsersicflux_error(data[i].DISKFIT,$
                                                                        data[i].COVAR_DISK,$
                                                                        /disk, cutoff=1)
                table[i].TOTAL_FLUX_ERR[0,color] = totalsersicflux_error(data[i].BULGEFIT,$
                                                                         data[i].COVAR_BULGE,$
                                                                         /bulge,/disk,cutoff=1)
                table[i].TOTAL_FLUX_ERR[1,color] = totalsersicflux_error(data[i].DISKFIT,$
                                                                         data[i].COVAR_DISK,$
                                                                         /bulge,/disk,cutoff=1)
                table[i].TOTAL_FLUX_ERR[4,color] = totalsersicflux_error(data[i].SERSICFIT,$
                                                                         data[i].COVAR_SERSIC,$
                                                                         cutoff=1)
                table[i].TOTAL_FLUX_ERR[2,color] = totalsersicflux_error(data[i].DVCFIT,$
                                                                         data[i].COVAR_DVC,$
                                                                         cutoff=1)
                table[i].TOTAL_FLUX_ERR[3,color] = totalsersicflux_error(data[i].EXPFIT,$
                                                                         data[i].COVAR_EXP,$
                                                                         cutoff=1)

                fpObj = mrdfits( sdss_name('fpObjc', data[i].RUN,$
                                           data[i].camcol, $
                                           data[i].field, $
                                           rerun=data[i].rerun), 1, $
                                 rows=[data[i].id-1], /silent)
                junk=mrdfits( sdss_name('fpObjc', data[i].RUN,$
                                        data[i].camcol, $
                                        data[i].field, $
                                        rerun=data[i].rerun), 0, $
                              header, /silent)
                table[i].SKY_COUNTS[*] = fpObj.SKY[*]
                table[i].SKY_ERR[*] = fpObj.SKYERR[*]
                table[i].GAIN = double((strsplit(fxpar(header, 'GAIN'), /extract)))


                ;get the central pixel and the min/max b-d
                makegrid, table[i].ATLAS_X, table[i].ATLAS_Y, x, y
                psffile = string( data[i].RERUN, data[i].RUN,$
                                  data[i].CAMCOL, data[i].RUN, $
                                  data[i].CAMCOL, data[i].FIELD, $
                                  format='("/u/dss/redux/",i0,"/",i0,"/objcs/",i0,"/psField-",i06,"-",i0,"-",i04,".fit")' )
                locpsf = sdss_psf_recon(mrdfits(psffile, color+1, /silent), $
                                        data[i].ROWC+0.5, $
                                        data[i].COLC+0.5, normalize=1.0 )


                bulge = pixelfluxpsf(x,y,[data[i].BULGEFIT[6:10],$
                                          data[i].BULGEFIT[12:13],$
                                          data[i].BULGEFIT[11]],$
                                     _EXTRA={cutoff:1,$
                                             psfImage:locpsf})
                disk =  pixelfluxpsf(x,y,[data[i].BULGEFIT[0:4],$
                                          data[i].BULGEFIT[12:13],$
                                          data[i].BULGEFIT[5]],$
                                     _EXTRA={cutoff:1,$
                                             psfImage:locpsf})
                table[i].BULGE_CENPIX[0] = max(bulge)
                table[i].DISK_CENPIX[0] = max(disk)
                table[i].B_MINUS_D_MAX[0] = max(bulge-disk,min=mn)
                table[i].B_MINUS_D_MIN[0] = mn

                bulge = pixelfluxpsf(x,y,[data[i].DISKFIT[6:10],$
                                          data[i].DISKFIT[12:13],$
                                          data[i].DISKFIT[11]],$
                                     _EXTRA={cutoff:1,$
                                             psfImage:locpsf})
                disk =  pixelfluxpsf(x,y,[data[i].DISKFIT[0:4],$
                                          data[i].DISKFIT[12:13],$
                                          data[i].DISKFIT[5]],$
                                     _EXTRA={cutoff:1,$
                                             psfImage:locpsf})
                table[i].BULGE_CENPIX[1] = max(bulge)
                table[i].DISK_CENPIX[1] = max(disk)
                table[i].B_MINUS_D_MAX[1] = max(bulge-disk,min=mn)
                table[i].B_MINUS_D_MIN[1] = mn

                bulge = pixelfluxpsf(x,y,data[i].DVCFIT,$
                                     _EXTRA={cutoff:1,$
                                             psfImage:locpsf})
                table[i].BULGE_CENPIX[2] = max(bulge)
                bulge = pixelfluxpsf(x,y,data[i].SERSICFIT,$
                                     _EXTRA={cutoff:1,$
                                             psfImage:locpsf})
                table[i].BULGE_CENPIX[4] = max(bulge)
                disk = pixelfluxpsf(x,y,data[i].EXPFIT,$
                                    _EXTRA={cutoff:1,$
                                            psfImage:locpsf})
                table[i].DISK_CENPIX[4] = max(disk)
                
                continue
                
            endif
                   
            atlas = get_imivarpsf(data[i].RUN, data[i].CAMCOL,$
                                  data[i].FIELD,data[i].ID,$
                                  data[i].ROWC[color],data[i].COLC[color],$
                                  rerun=data[i].RERUN,filter=color)
                      
            table[i].ATLAS_X = (size(atlas.image, /dimensions))[0]
            table[i].ATLAS_Y = (size(atlas.image, /dimensions))[1]
            if not is_struct(atlas) then begin
                print, 'no atlas image for ',data[i].RUN,data[i].CAMCOL,$
                  data[i].FIELD,data[i].ID,data[i].RERUN
                continue
            endif

            jj = where(atlas.ivar lt 1.e-10, nn_ivar)
            base_dof = n_elements(atlas.image) - nn_ivar 

            bb = scalefit(atlas.image,atlas.psf,atlas.ivar, $
                          switch_band(data[i].BULGEFIT[*,rband],2,color,data[i]), covar=bbcovar)

            ;table[i].CHI_SQ_SCALED[0,color]  = getChiSquare(atlas.image, atlas.ivar, bb, $
            ;                                                _EXTRA={psfImage:atlas.psf, cutoff:1})/$
            ;  (base_dof - 2.)
            table[i].DOF_SCALED[0,color] = base_dof-2
            
            if(data[i].BULGEFIT[6,rband] ne 0.0 ) then table[i].BULGE_REL_RBAND[0,color] = bb[6]/ $
              data[i].BULGEFIT[6,rband] $
            else table[i].BULGE_REL_RBAND[0,color] = 0.0
            if(data[i].BULGEFIT[0,rband] ne 0.0 ) then table[i].DISK_REL_RBAND[0,color] = bb[0]/$
              data[i].BULGEFIT[0,rband] $
            else table[i].DISK_REL_RBAND[0,color] = 0.0
    
            dd = scalefit(atlas.image,atlas.psf,atlas.ivar,$
                          switch_band(data[i].DISKFIT[*,rband],2,color,data[i]), covar=ddcovar)
            ;table[i].CHI_SQ_SCALED[1,color]  = getChiSquare(atlas.image, atlas.ivar, dd, $
            ;                                                _EXTRA={psfImage:atlas.psf, cutoff:1})/ $
            ;  (base_dof - 2.)
            table[i].DOF_SCALED[1,color] = base_dof-2
            
            if(data[i].DISKFIT[6,rband] ne 0.0 ) then table[i].BULGE_REL_RBAND[1,color] = dd[6] / $
              data[i].DISKFIT[6,rband] $
            else table[i].BULGE_REL_RBAND[1,color] = 0.0
            if(data[i].DISKFIT[0,rband] ne 0.0 ) then table[i].DISK_REL_RBAND[1,color] = dd[0] / $
              data[i].DISKFIT[0,rband] $
            else table[i].DISK_REL_RBAND[1,color] = 0.0

            dvc = scalefit(atlas.image,atlas.psf,atlas.ivar, $
                           switch_band(data[i].DVCFIT[*,rband], 2,color,data[i]),covar=dvccovar)
            exp = scalefit(atlas.image,atlas.psf,atlas.ivar, $
                           switch_band(data[i].EXPFIT[*,rband],2,color,data[i]),covar=expcovar)
            ser = scalefit(atlas.image,atlas.psf,atlas.ivar, $
                           switch_band(data[i].SERSICFIT[*,rband],2,color,data[i]),covar=sercovar)
            ;table[i].CHI_SQ_SCALED[2,color]  = getChiSquare(atlas.image, atlas.ivar, dvc, $
            ;                                                _EXTRA={psfImage:atlas.psf, cutoff:1}) / $
            ;  (base_dof - 1.)
            table[i].DOF_SCALED[2,color] = base_dof-1
            ;table[i].CHI_SQ_SCALED[3,color]  = getChiSquare(atlas.image, atlas.ivar, exp, $
            ;                                                _EXTRA={psfImage:atlas.psf, cutoff:1}) / $
            ;  (base_dof - 1.)
            table[i].DOF_SCALED[3,color] = base_dof-1
            ;table[i].CHI_SQ_SCALED[4,color]  = getChiSquare(atlas.image, atlas.ivar, ser, $
            ;                                                _EXTRA={psfImage:atlas.psf, cutoff:1}) / $
            ;  (base_dof - 1.)
            table[i].DOF_SCALED[4,color] = base_dof-1
            
           ;get rid of header pointers
            ;heap_free, phdr

            if( data[i].DVCFIT[0,rband] ne 0.0) then table[i].BULGE_REL_RBAND[2,color] = dvc[0]/$
              data[i].DVCFIT[0,rband] $
            else table[i].BULGE_REL_RBAND[2,color] = 0.0
            if(data[i].EXPFIT[0,rband] ne 0.0) then table[i].DISK_REL_RBAND[3,color] = exp[0]/$
              data[i].EXPFIT[0,rband] $
            else table[i].DISK_REL_RBAND[3,color] = 0.0
            if( data[i].SERSICFIT[0,rband] ne 0.0) then table[i].BULGE_REL_RBAND[4,color] = ser[0]/$
              data[i].SERSICFIT[0,rband] $
            else table[i].BULGE_REL_RBAND[4,color] = 0.0

            table[i].BULGE_FLUX_ERR[0,color]  = get_color_error(data[i].BULGEFIT,$
                                                                data[i].COVAR_BULGE,$
                                                                [table[i].DISK_REL_RBAND[0,color],$
                                                                 table[i].BULGE_REL_RBAND[0,color]],$
                                                                bbcovar,/bulge, cutoff=1)
            table[i].BULGE_FLUX_ERR[1,color] =  get_color_error(data[i].DISKFIT,$
                                                                data[i].COVAR_DISK,$
                                                                [table[i].DISK_REL_RBAND[1,color],$
                                                                 table[i].BULGE_REL_RBAND[1,color]],$
                                                                ddcovar,/bulge, cutoff=1)
            table[i].DISK_FLUX_ERR[0,color]  = get_color_error(data[i].BULGEFIT,$
                                                               data[i].COVAR_BULGE,$
                                                               [table[i].DISK_REL_RBAND[0,color],$
                                                                table[i].BULGE_REL_RBAND[0,color]],$
                                                               bbcovar,/disk, cutoff=1)
            table[i].DISK_FLUX_ERR[1,color] =  get_color_error(data[i].DISKFIT,$
                                                               data[i].COVAR_DISK,$
                                                               [table[i].DISK_REL_RBAND[1,color],$
                                                                 table[i].BULGE_REL_RBAND[1,color]],$
                                                               ddcovar,/disk, cutoff=1)
            table[i].TOTAL_FLUX_ERR[0,color] = get_color_error(data[i].BULGEFIT,$
                                                                data[i].COVAR_BULGE,$
                                                                [table[i].DISK_REL_RBAND[0,color],$
                                                                 table[i].BULGE_REL_RBAND[0,color]],$
                                                                bbcovar, cutoff=1)
            table[i].TOTAL_FLUX_ERR[1,color] = get_color_error(data[i].DISKFIT,$
                                                               data[i].COVAR_DISK,$
                                                               [table[i].DISK_REL_RBAND[1,color],$
                                                                table[i].BULGE_REL_RBAND[1,color]],$
                                                               ddcovar, cutoff=1)
            table[i].TOTAL_FLUX_ERR[4,color] = get_color_error(data[i].SERSICFIT,$
                                                               data[i].COVAR_SERSIC,$
                                                               table[i].BULGE_REL_RBAND[4,color],$
                                                               sercovar, cutoff=1)
            table[i].TOTAL_FLUX_ERR[2,color] = get_color_error(data[i].DVCFIT,$
                                                               data[i].COVAR_DVC,$
                                                               table[i].BULGE_REL_RBAND[2,color],$
                                                               dvccovar, cutoff=1)
            table[i].TOTAL_FLUX_ERR[3,color] = get_color_error(data[i].EXPFIT,$
                                                               data[i].COVAR_EXP,$
                                                               table[i].DISK_REL_RBAND[3,color],$
                                                               expcovar, cutoff=1)
            table[i].BULGE_FLUX_ERR[2,color] = table[i].TOTAL_FLUX_ERR[2,color]
            table[i].BULGE_FLUX_ERR[4,color] = table[i].TOTAL_FLUX_ERR[4,color]
            table[i].DISK_FLUX_ERR[3,color] = table[i].TOTAL_FLUX_ERR[3,color]
          
        endfor 
    endif else begin
        table[i].BULGE_REL_RBAND=0.0
        table[i].DISK_REL_RBAND=0.0
        my_color = (where(bands eq fiducial))[0]
        table[i].BULGE_REL_RBAND[*,my_color]=1.
        table[i].DISK_REL_RBAND[*,my_color]=1.
        table[i].BULGE_FLUX_ERR[0,my_color]  = totalsersicflux_error(data[i].BULGEFIT,$
                                                                  data[i].COVAR_BULGE,$
                                                                  /bulge, cutoff=1)
        table[i].BULGE_FLUX_ERR[1,my_color] = totalsersicflux_error(data[i].DISKFIT,$
                                                                 data[i].COVAR_DISK,$
                                                                 /bulge, cutoff=1)
        table[i].DISK_FLUX_ERR[0,my_color]  = totalsersicflux_error(data[i].BULGEFIT,$
                                                                 data[i].COVAR_BULGE,$
                                                                 /disk, cutoff=1)
        table[i].DISK_FLUX_ERR[1,my_color] = totalsersicflux_error(data[i].DISKFIT,$
                                                                data[i].COVAR_DISK,$
                                                                /disk, cutoff=1)
        table[i].TOTAL_FLUX_ERR[0,my_color] = totalsersicflux_error(data[i].BULGEFIT,$
                                                                 data[i].COVAR_BULGE,$
                                                                 /bulge,/disk,cutoff=1)
        table[i].TOTAL_FLUX_ERR[1,my_color] = totalsersicflux_error(data[i].DISKFIT,$
                                                                 data[i].COVAR_DISK,$
                                                                 /bulge,/disk,cutoff=1)
        table[i].TOTAL_FLUX_ERR[4,my_color] = totalsersicflux_error(data[i].SERSICFIT,$
                                                                 data[i].COVAR_SERSIC,$
                                                                 cutoff=1)
        table[i].TOTAL_FLUX_ERR[2,my_color] = totalsersicflux_error(data[i].DVCFIT,$
                                                                 data[i].COVAR_DVC,$
                                                                 cutoff=1)
        table[i].TOTAL_FLUX_ERR[3,my_color] = totalsersicflux_error(data[i].EXPFIT,$
                                                                    data[i].COVAR_EXP,$
                                                                    cutoff=1)
        table[i].BULGE_FLUX_ERR[2,my_color] = table[i].TOTAL_FLUX_ERR[2,my_color]
        table[i].BULGE_FLUX_ERR[4,my_color] = table[i].TOTAL_FLUX_ERR[4,my_color]
        table[i].DISK_FLUX_ERR[3,my_color] = table[i].TOTAL_FLUX_ERR[3,my_color]

                                                                                 
        fpObj = mrdfits( sdss_name('fpObjc', data[i].RUN,$
                                   data[i].camcol, $
                                   data[i].field, $
                                   rerun=data[i].rerun), 1, $
                         rows=[data[i].id-1], /silent)
        junk=mrdfits( sdss_name('fpObjc', data[i].RUN,$
                                data[i].camcol, $
                                data[i].field, $
                                rerun=data[i].rerun), 0, $
                      header, /silent)
        table[i].SKY_COUNTS[*] = fpObj.SKY[*]
        table[i].SKY_ERR[*] = fpObj.SKYERR[*]
        table[i].GAIN = double((strsplit(fxpar(header, 'GAIN'), /extract)))
               
        atlas = get_imivarpsf(data[i].RUN, data[i].CAMCOL,$
                              data[i].FIELD,data[i].ID,$
                              data[i].ROWC[my_color],data[i].COLC[my_color],$
                              rerun=data[i].RERUN,filter=my_color)

        table[i].ATLAS_X = (size(atlas.image, /dimensions))[0]
        table[i].ATLAS_Y = (size(atlas.image, /dimensions))[1]
        makegrid, table[i].ATLAS_X, table[i].ATLAS_Y, x, y
        
        bulge = pixelfluxpsf(x,y,[data[i].BULGEFIT[6:10],$
                                  data[i].BULGEFIT[12:13],$
                                  data[i].BULGEFIT[11]],$
                             _EXTRA={cutoff:1,$
                                     psfImage:atlas.psf})
        disk =  pixelfluxpsf(x,y,[data[i].BULGEFIT[0:4],$
                                  data[i].BULGEFIT[12:13],$
                                  data[i].BULGEFIT[5]],$
                             _EXTRA={cutoff:1,$
                                     psfImage:atlas.psf})
        table[i].BULGE_CENPIX[0] = max(bulge)
        table[i].DISK_CENPIX[0] = max(disk)
        table[i].B_MINUS_D_MAX[0] = max(bulge-disk,min=mn)
        table[i].B_MINUS_D_MIN[0] = mn
        
        bulge = pixelfluxpsf(x,y,[data[i].DISKFIT[6:10],$
                                  data[i].DISKFIT[12:13],$
                                  data[i].DISKFIT[11]],$
                             _EXTRA={cutoff:1,$
                                     psfImage:atlas.psf})
        disk =  pixelfluxpsf(x,y,[data[i].DISKFIT[0:4],$
                                  data[i].DISKFIT[12:13],$
                                  data[i].DISKFIT[5]],$
                             _EXTRA={cutoff:1,$
                                     psfImage:atlas.psf})
        table[i].BULGE_CENPIX[1] = max(bulge)
        table[i].DISK_CENPIX[1] = max(disk)
        table[i].B_MINUS_D_MAX[1] = max(bulge-disk,min=mn)
        table[i].B_MINUS_D_MIN[1] = mn
        
        bulge = pixelfluxpsf(x,y,data[i].DVCFIT,$
                             _EXTRA={cutoff:1,$
                                     psfImage:atlas.psf})
        table[i].BULGE_CENPIX[2] = max(bulge)
        bulge = pixelfluxpsf(x,y,data[i].SERSICFIT,$
                             _EXTRA={cutoff:1,$
                                     psfImage:atlas.psf})
        table[i].BULGE_CENPIX[4] = max(bulge)
        disk = pixelfluxpsf(x,y,data[i].EXPFIT,$
                            _EXTRA={cutoff:1,$
                                    psfImage:atlas.psf})
        table[i].DISK_CENPIX[4] = max(disk)

    endelse

    ;r-band flux
    rn4BulgeParams = [data[i].BULGEFIT[6:10,rband],data[i].BULGEFIT[12:13,rband],$
                      data[i].BULGEFIT[11,rband]]
    rn4DiskParams = [data[i].BULGEFIT[0:4,rband],data[i].BULGEFIT[12:13,rband],$
                     data[i].BULGEFIT[5,rband]]
    if( rn4BulgeParams[3] ne -1.0 ) then $
      rn4Bulge = totalsersicflux(rn4BulgeParams,cutoff=1) $
    else rn4Bulge = -1.0
    if( rn4DiskParams[3] ne -1.0 ) then $
      rn4Disk = totalsersicflux(rn4DiskParams,cutoff=1) $
    else rn4Disk = -1.0

    ;fluxes in Maggies
    f1 = table[i].NMGYPERCOUNT[*] * 1.0e-9 * $
      table[i].BULGE_REL_RBAND[0,*] * rn4Bulge
    f2 = table[i].NMGYPERCOUNT[*] * 1.0e-9 * $
      table[i].DISK_REL_RBAND[0,*] * rn4Disk
    
    table[i].BULGE_FLUX_COUNTS[0,*] = f1/(1.e-9*table[i].NMGYPERCOUNT[*])
    table[i].DISK_FLUX_COUNTS[0,*] = f2/(1.e-9*table[i].NMGYPERCOUNT[*])
    table[i].TOTAL_FLUX_COUNTS[0,*] = (f1+f2)/(1.e-9*table[i].NMGYPERCOUNT[*])
    
    table[i].BULGE_MAG[0,*] = sdss_asinh_mag(f1,filter=[0,1,2,3,4])
    table[i].DISK_MAG[0,*] = sdss_asinh_mag(f2,filter=[0,1,2,3,4])
    table[i].TOTAL_MAG[0,*] = sdss_asinh_mag(f1+f2,filter=[0,1,2,3,4])
    
    rn1BulgeParams = [data[i].DISKFIT[6:10,rband],data[i].DISKFIT[12:13,rband],$
                      data[i].DISKFIT[11,rband]]
    rn1DiskParams = [data[i].DISKFIT[0:4,rband],data[i].DISKFIT[12:13,rband],$
                     data[i].DISKFIT[5,rband]]
    if( rn1BulgeParams[3] ne -1.0 ) then $
      rn1Bulge=totalsersicflux(rn1BulgeParams,cutoff=1) $
    else rn1Bulge = -1.0
    if( rn1DiskParams[3] ne -1.0 ) then $
      rn1Disk = totalsersicflux(rn1DiskParams,cutoff=1) $
    else rn1Disk = -1.0
    ;fluxes in Maggies
    f1 = table[i].NMGYPERCOUNT[*] * 1.0e-9 * $
      table[i].BULGE_REL_RBAND[1,*] * rn1Bulge
    f2 = table[i].NMGYPERCOUNT[*] * 1.0e-9 * $
      table[i].DISK_REL_RBAND[1,*] * rn1Disk

    table[i].BULGE_FLUX_COUNTS[1,*] = f1/(1.e-9*table[i].NMGYPERCOUNT[*])
    table[i].DISK_FLUX_COUNTS[1,*] = f2/(1.e-9*table[i].NMGYPERCOUNT[*])
    table[i].TOTAL_FLUX_COUNTS[1,*] = (f1+f2)/(1.e-9*table[i].NMGYPERCOUNT[*])
 
    table[i].BULGE_MAG[1,*] = sdss_asinh_mag(f1,filter=[0,1,2,3,4])
    table[i].DISK_MAG[1,*] = sdss_asinh_mag(f2,filter=[0,1,2,3,4])
    table[i].TOTAL_MAG[1,*] = sdss_asinh_mag(f1+f2,filter=[0,1,2,3,4])

    dvcbulge = totalsersicflux(data[i].DVCFIT[*,rband],cutoff=1)
    expdisk = totalsersicflux(data[i].EXPFIT[*,rband],cutoff=1)
    serbulge = totalsersicflux(data[i].SERSICFIT[*,rband],cutoff=1)
    ;fluxes in Maggies
    f1 = table[i].NMGYPERCOUNT[*] * 1.0e-9 * $
      table[i].BULGE_REL_RBAND[2,*] * dvcbulge
    f2 = table[i].NMGYPERCOUNT[*] * 1.0e-9 * $
      table[i].DISK_REL_RBAND[3,*] * expdisk
    f3 = table[i].NMGYPERCOUNT[*] * 1.0e-9 * $
      table[i].BULGE_REL_RBAND[4,*] * serbulge

    table[i].BULGE_FLUX_COUNTS[2,*] = f1/(1.e-9*table[i].NMGYPERCOUNT[*])
    table[i].TOTAL_FLUX_COUNTS[2,*] = (f1)/(1.e-9*table[i].NMGYPERCOUNT[*])
    table[i].DISK_FLUX_COUNTS[3,*] = f2/(1.e-9*table[i].NMGYPERCOUNT[*])
    table[i].TOTAL_FLUX_COUNTS[3,*] = (f2)/(1.e-9*table[i].NMGYPERCOUNT[*])
    table[i].BULGE_FLUX_COUNTS[4,*] = f3/(1.e-9*table[i].NMGYPERCOUNT[*])
    table[i].TOTAL_FLUX_COUNTS[4,*] = (f3)/(1.e-9*table[i].NMGYPERCOUNT[*])

    table[i].BULGE_MAG[2,*] = sdss_asinh_mag(f1,filter=[0,1,2,3,4])
    table[i].TOTAL_MAG[2,*] = sdss_asinh_mag(f1,filter=[0,1,2,3,4])
    table[i].DISK_MAG[3,*] = sdss_asinh_mag(f2,filter=[0,1,2,3,4])
    table[i].TOTAL_MAG[3,*] = sdss_asinh_mag(f2,filter=[0,1,2,3,4])
    table[i].BULGE_MAG[4,*] = sdss_asinh_mag(f3,filter=[0,1,2,3,4])
    table[i].TOTAL_MAG[4,*] = sdss_asinh_mag(f3,filter=[0,1,2,3,4])



    
    for f=0, n_elements(fitbands)-1 do begin
        filter = data[0].BANDS[f]
        bulgeToTotn4=bulgetotot(data[i].BULGEFIT[*,f], cutoff=1)
        
        bulgeToTotn1 = bulgetotot(data[i].DISKFIT[*,f], cutoff=1)
        
        table[i].BULGE_TO_TOT[0,fitbands[f]] = bulgeToTotn4
        table[i].BULGE_TO_TOT[1,fitbands[f]] = bulgeToTotn1
                                        ;bulge to total errors:
        table[i].BULGE_TO_TOT_ERR[0,fitbands[f]] = $
          bulgetotot_error( data[i].BULGEFIT[*,f], $
                            data[i].COVAR_BULGE[*,*,f] )
        table[i].BULGE_TO_TOT_ERR[1,fitbands[f]] = $
          bulgetotot_error( data[i].DISKFIT[*,f], $
                            data[i].COVAR_DISK[*,*,f] )

        

        
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                        ;make the output images
        if( keyword_set(images) ) then begin
            atlas = (sdss_atlas_image( data[i].RUN, $
                                       data[i].CAMCOL, $
                                       data[i].FIELD, $
                                       data[i].ID, $
                                       rerun=data[i].RERUN, $
                                       filter=filter ))*1.0D
            if n_elements(atlas le 1) then begin
                print, 'no atlas image for ',data[i].RUN, data[i].CAMCOL,$
                  data[i].FIELD, data[i].ID, data[i].RERUN
                continue
            endif
            ;path="../z_rescale/S0013-00941_zmove/"+data[i].GALID+"-Z-"
            ;red = where(abs(data[i].MOVED_Z-data[i].CURRENT_Z) lt 1.0e-5)
            ;print, red
            ;im_file=path+string(red,format='(i02,".fits")')

            ;atlas = mrdfits(im_file,1,/silent)
            psffile = sdss_name( 'psField', data[i].RUN, $
                                 data[i].CAMCOL, data[i].FIELD, $
                                 rerun=data[i].RERUN )
            locpsf = sdss_psf_recon(mrdfits(psffile, filter+1, /silent), $
                                    data[i].ROWC[filter]+0.5, $
                                    data[i].COLC[filter]+0.5, normalize=1.0 )
            
            lenx = (size(atlas, /dimensions))[0]
            leny = (size(atlas, /dimensions))[1]
            
            makegrid, lenx, leny, x, y
            
            bulgedisk = pixelfluxpsf( x, y, data[i].BULGEFIT[*,f], $
                                      _EXTRA={psfImage:locpsf,cutoff:1})+1000.0
            bulge = pixelfluxpsf(x,y,[data[i].BULGEFIT[6:10], $
                                      data[i].BULGEFIT[12:13], $
                                      data[i].BULGEFIT[11]], $
                                 _EXTRA={psfImage:locpsf,cutoff:1})
                              
            bd_nopsf = pixelfluxpsf(x,y,data[i].BULGEFIT[*,f], $
                                    _EXTRA={cutoff:1})+1000.0
            diskdisk = pixelfluxpsf( x, y, data[i].DISKFIT[*,f], $
                                     _EXTRA={psfImage:locpsf,cutoff:1})+1000.0
            dd_nopsf = pixelfluxpsf(x,y,data[i].DISKFIT[*,f], $
                                    _EXTRA={cutoff:1})+1000.0
            disk = pixelfluxpsf(x,y,[data[i].DISKFIT[6:10], $
                                     data[i].DISKFIT[12:13], $
                                     data[i].DISKFIT[11]], $
                                _EXTRA={psfImage:locpsf,cutoff:1})
            imagesn4 = inputfolder+'/images/'+table[i].GALID+'-n4.fits'
;               string(data[i].CURRENT_Z,format='(F5.3)')+'-n4.fits'
            mwrfits, bulgedisk, imagesn4, ['IMAGE = model_'+strtrim(filter,2),''], /create
            mwrfits, bulgedisk-atlas, imagesn4, ['IMAGE = model-atlas_'+strtrim(filter,2),'']
            mwrfits, bd_nopsf, imagesn4, ['IMAGE = model_nopsf_'+strtrim(filter,2),'']
            mwrfits, bulge, imagesn4, ['IMAGE = bulgemodel'+strtrim(filter,2),'']
            imagesn1=inputfolder+'/images/'+table[i].GALID+'-n1.fits'
 ;             string(data[i].CURRENT_Z,format='(F5.3)')+'-n1.fits'
            mwrfits, diskdisk, imagesn1, ['IMAGE = model_'+strtrim(filter,2),''], /create
            mwrfits, diskdisk-atlas, imagesn1, ['IMAGE = model-atlas_'+strtrim(filter,2),'']
            mwrfits, dd_nopsf, imagesn1, ['IMAGE = model_nopsf+'+strtrim(filter,2),'']
            mwrfits, disk, imagesn1, ['IMAGE = bulgemodel'+strtrim(filter,2),'']
            ;one component fits
            imagesexp=inputfolder+'/images/'+table[i].GALID+'-exp.fits'
            disk = pixelfluxpsf(x,y,data[i].EXPFIT[*,f], $
                                _EXTRA={psfImage:locpsf,cutoff:1})+1000.
            mwrfits, disk, imagesexp, ['IMAGE = model_'+strtrim(filter,2),''], /create
            mwrfits, disk-atlas, imagesexp, ['IMAGE = model-atlas_'+strtrim(filter,2),'']
            imagesdvc=inputfolder+'/images/'+table[i].GALID+'-dvc.fits'
            disk = pixelfluxpsf(x,y,data[i].DVCFIT[*,f], $
                                _EXTRA={psfImage:locpsf,cutoff:1})+1000.
            mwrfits, disk, imagesdvc, ['IMAGE = model_'+strtrim(filter,2),''], /create
            mwrfits, disk-atlas, imagesdvc, ['IMAGE = model-atlas_'+strtrim(filter,2),'']
            imagesser=inputfolder+'/images/'+table[i].GALID+'-ser.fits'
            disk = pixelfluxpsf(x,y,data[i].SERSICFIT[*,f], $
                                _EXTRA={psfImage:locpsf,cutoff:1})+1000.
            mwrfits, disk, imagesser, ['IMAGE = model_'+strtrim(filter,2),''], /create
            mwrfits, disk-atlas, imagesser, ['IMAGE = model-atlas_'+strtrim(filter,2),'']

            ;sdss softening parameters
            soft=[1.4e-10,0.9e-10,1.2e-10,1.8e-10,7.4e-10]
            ps_plot, atlas, bulgedisk, diskdisk, $
              table[i].GALID,$;+string(data[i].CURRENT_Z,format='(F5.3)'), $
              bands[filter], table[i].NMGYPERCOUNT[filter], soft[filter], inputfolder
        endif
        
    
    endfor
endfor

mwrfits, table, outputtable, /create

;do the kcorrection here
add_kcorrect, outputtable, outputRAW

END
