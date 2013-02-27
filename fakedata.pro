;____________________________________
;
; Claire Lackner
;  Jan 20 2010
;
;  fakedata.pro
;
;  makes a suite of fake galaxy images and fits at a variety of radii
;   writing the final info to a fits file
; 
;____________________________________



PRO fakedata, numGal, outfile

seed = systime(1)
prof = floor( randomu(seed, numGal)*3.0+1.0 )
fit_range = [5.0, 1.0, 0.6, 0.5, 0.3, 0.2, 0.1, 0.05]
nfits = 9
outputentry = { PROFILE:prof[0], FITRANGE:dblarr(9), NORMS:dblarr(11), SERSN_R0:dblarr(11), $
                R50_H0:dblarr(11), ELLIP_INC:dblarr(11), C_RMAX:dblarr(11), $
                XCEN:dblarr(11), YCEN:dblarr(11), PHI:dblarr(11), $
                CHI_SQUARE:dblarr(11), PSF_FWHM:0.0 }
output = replicate( outputentry, numGal )

val = 0
psf = 0
image = 0

for i=0, numGal-1 do begin
    fitfakes, prof[i], val, fit_range, image

    output[i].PROFILE = prof[i]
    output[i].FITRANGE = fit_range
    output[i].NORMS = val[*,0]
    output[i].SERSN_R0 = val[*,1]
    output[i].R50_H0 = val[*,2]
    output[i].ELLIP_INC = val[*,3]
    output[i].C_RMAX = val[*,4]
    output[i].XCEN = val[*,5]
    output[i].YCEN = val[*,6]
    output[i].PHI = val[*,7]
    output[i].CHI_SQUARE=val[*,9]
    output[i].PSF_FWHM = val[0,8]

    ;write image to a file
    outname=string(i, prof[i], format='("fake",i04,"-",i1,".fits")')
    mwrfits, image, outname
  
endfor


mwrfits, output, outfile


END


