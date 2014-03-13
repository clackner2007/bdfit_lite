;+
;;;;;
;Claire Lackner
;March 13 2014
;
; get_sdsspsf.pro
;
; function to return the sdss psf for a given run, rerun, camcol,
;field, rowc, colc, and filter
; filter is given as 'u', 'g', 'r', 'i', or 'z'
;
; 
; returns as a normalized array
;
;-

FUNCTION get_sdsspsf, run, rerun, cacmol, field, rowc, colc, filter=filter

  bands = ['u', 'g', 'r', 'i', 'z']
  band = bands[filter]

  psffile = string( rerun, run, camcol, run, camcol, field, $
                    format='("/u/dss/redux/",i0,"/",i0,"/objcs/",i0,"/psField-",i06,"-",i0,"-",i04,".fit")' )

  locpsf = sdss_psf_recon(mrdfits(psffile, band+1, /silent), $
                         rowc+0.5, colc+0.5, normalize=1.0)


END
