;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; switch_band.pro
;
; Claire Lackner
;
;
; FUNCTION SWITCH_BAND, params, band_in, band_out
;
; RETURNS the parameters for an atlas image fit in a different band
;  changes the x-y coordinates of the center using align_atlas_point
;  modifies the the rotation angles by phi_offset
;
;
;  needs an sdss_object containing: RUN, CAMCOL, FIELD, ID, RERUN, PHI_OFFSET(radians)
; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION switch_band, params, band_in, band_out, sdss_object, $
                      fpc_file=fpc_file, fpcheader=fpcheader

new_params = params

if band_in eq band_out then return, new_params

if( n_elements(new_params) eq 8 ) then begin
    x_p=5
    y_p=6
    phi_p=[7]
endif else begin
    x_p=12
    y_p=13
    phi_p=[5,11]
endelse

;align_atlas_point, params[x_p], params[y_p], xx, yy, band_in, band_out, $
;  sdss_object.RUN, sdss_object.CAMCOL, sdss_object.FIELD, sdss_object.ID, $
;  rerun=sdss_object.RERUN, fpc_file=fpc_file, fpcheader=fpcheader
align_atlas_astrans, params[x_p], params[y_p], xx, yy, band_in, band_out, $
  sdss_object.RUN, sdss_object.CAMCOL, sdss_object.FIELD, sdss_object.ID, $
  rerun=sdss_object.RERUN


new_params[x_p] = xx
new_params[y_p] = yy

;rotate for new phi offset
new_params[phi_p] = params[phi_p] - sdss_object.PHI_OFFSET_RAD[band_in] + $
  sdss_object.PHI_OFFSET_RAD[band_out]


return, new_params

END
