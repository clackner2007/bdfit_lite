;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
;
;  makes postscript plots of image,fit and residuals
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION getluptitude, counts, nmgypercount, soft

return, $
  (asinh(counts*1.0e-9*nmgypercount/(2.*soft)) + alog(soft))*(2.5/alog(10.0))

END


PRO ps_plot, atlas, fitn4, fitn1, galid, band, nmgypercount, softening, folder

set_plot, "ps"
!P.MULTI = [0,5,1]
!P.CHARSIZE=1.8
output=folder+'/images/'+galid+'-'+band+'.eps'
device, /encapsulated, filename=output, /cm, landscape=0
device, bits_per_pixel=8
device, xsize= 30.0

!P.TITLE = galid+'-'+band
;tvscale, atlas, /keep_aspect_ratio, multimargin=[1,0,2,0]
tvscale, getluptitude(atlas-1000.,nmgypercount, softening),$
  /keep_aspect_ratio, multimargin=[1,0,2,0]
!P.TITLE='n=4'
tvscale, getluptitude(fitn4-1000.,nmgypercount, softening), $
  /keep_aspect_ratio, multimargin=[1,0,2,0]
!P.TITLE='n=4,resid'
tvscale, getluptitude(fitn4-atlas,nmgypercount, softening), $
  /keep_aspect_ratio, multimargin=[1,0,2,0]
!P.TITLE='n=1'
tvscale, getluptitude(fitn1-1000.,nmgypercount, softening),$
  /keep_aspect_ratio, multimargin=[1,0,2,0]
!P.TITLE='n=1,resid'
tvscale, getluptitude(fitn1-atlas,nmgypercount, softening), $
  /keep_aspect_ratio, multimargin=[1,0,2,0]

device, /close
set_plot, "x"
!P.TITLE=''
!P.MULTI=0
!Y.OMARGIN = [0,0]
!P.CHARSIZE=1.0
END
