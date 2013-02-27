;;;;;;;;;;;;;;;;;;;;
;;; returns the error in the scaled color
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION get_color_error, fid_fit, covar_fid, scaling, scaling_covar, $
                          cutoff=cutoff,bulge=bulge, disk=disk

;return error for the bulge
if keyword_set(bulge) then begin
    bulge_orig = totalsersicflux([fid_fit[6:10],0.,0.,0.,0.],cutoff=cutoff)
    bulge_orig_err = totalsersicflux_error(fid_fit, covar_fid, /bulge, cutoff=cutoff)

    return, sqrt(bulge_orig^2*scaling_covar[3]+scaling[1]^2*bulge_orig_err^2)
endif else begin
    ;return error for the disk
    if keyword_set(disk) then begin
        disk_orig = totalsersicflux([fid_fit[0:4],0.,0.,0.,0.],cutoff=cutoff)
        disk_orig_err = totalsersicflux_error(fid_fit, covar_fid, /disk, cutoff=cutoff)
        
        return, sqrt(disk_orig^2*scaling_covar[0]+scaling[0]^2*disk_orig_err^2)
        ;return error for the total flux
    endif else begin
        tot_orig = totalsersicflux(fid_fit, cutoff=cutoff)
        
        if n_elements(fid_fit) eq 8 then begin
            tot_orig_error = totalsersicflux_error(fid_fit, covar_fid, cutoff=cutoff)
            return, sqrt( tot_orig^2*scaling_covar + scaling^2*tot_orig_error^2)
        endif else begin
            tot_orig_error = totalsersicflux_error(fid_fit, covar_fid, /disk, /bulge, cutoff=cutoff)
            bulge_orig = totalsersicflux([fid_fit[6:10],0.,0.,0.,0.],cutoff=cutoff)
            bulge_orig_err = totalsersicflux_error(fid_fit, covar_fid, /bulge, cutoff=cutoff)
            disk_orig = totalsersicflux([fid_fit[0:4],0.,0.,0.,0.],cutoff=cutoff)
            disk_orig_err = totalsersicflux_error(fid_fit, covar_fid, /disk, cutoff=cutoff)

            bd_covar = (tot_orig_error^2 - bulge_orig_err^2 - disk_orig_err^2) / $
              (2.*bulge_orig*disk_orig)

            return, sqrt(scaling_covar[3]*bulge_orig^2 + $
                        bulge_orig_err^2*scaling[1]^2 + $
                        scaling_covar[0]*disk_orig^2 + $
                        scaling[0]^2*disk_orig_err^2 + $
                        bd_covar*scaling[0]*scaling[1] + $
                        scaling_covar[1]*bulge_orig*disk_orig)
            
            
        endelse
    endelse
endelse 



END
