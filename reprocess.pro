;Claire Lackner
; Jul 16 2010
;
; reprocess.pro
;
; program to reprocess reina's fits to include the value for R_80
; containing 80% of the flux of the model
;

PRO reprocess, file, outputfile

data=mrdfits(file,1)
nfits = (size(data.BULGE_SERS_N))[1]
nbands = (size(data.BULGE_SERS_N))[2]

output_entry=create_struct( data[0], ['R_80'], dblarr(nfits,nbands) )
output = replicate(output_entry, n_elements(data) )
struct_assign, data, output




for gal=0,n_elements(data)-1 do begin
    params = [[[data[gal].DISK_IO/exp(kappa_sersic(1.0))]],$
          [[data[gal].DISK_SCALE/0.59684]],$
          [[dblarr(nfits,nbands) + 1.0]],$
          [[cos(data[gal].DISK_INC)]],$
          [[dblarr(nfits,nbands)]],$
          [[data[gal].DISK_PHI]],$
          [[data[gal].BULGE_IO/exp(kappa_sersic(data[gal].BULGE_SERS_N))]],$
          [[data[gal].BULGE_REFF]],$
          [[data[gal].BULGE_SERS_N]],$
          [[data[gal].BULGE_Q]],$
          [[dblarr(nfits,nbands)]],$
          [[data[gal].BULGE_PHI]],$
          [[data[gal].CEN_X]],[[data[gal].CEN_Y]]]
    for fit=0,nfits-1 do begin
        for band=0,nbands-1 do begin
            if( params[fit,band,13] eq 0 ) then r_80=0.0 $
            else begin
              if params[fit,band,8] le 1.0e-3 then $
                r_80 = get_r80([reform(params[fit,band,0:4]),0.0,0.0,params[fit,band,5]]) $
              else if params[fit,band,1] le 1.0e-3 then $
                r_80 = get_r80([reform(params[fit,band,6:10]),0.0,0.0,params[fit,band,11]]) $
              else r_80 = get_r80(reform(params[fit,band,*]))
          endelse
            output[gal].R_80[fit,band] = r_80
        endfor
    endfor
endfor

mwrfits, output, outputfile, /create

END
