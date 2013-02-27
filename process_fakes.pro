;idl code to add fluxes and errors to totalRAW


pro process_fakes, file_in, file_out

d=mrdfits(file_in,1)

oe = create_struct(d[0], ['TOTAL_FLUX', 'TOTAL_FLUX_ERR', 'BTT', 'BTT_ERR',$
                          'BULGE_FLUX', 'BULGE_FLUX_ERR', 'DISK_FLUX', "DISK_FLUX_ERR"],$
                   dblarr(5), dblarr(5), dblarr(5), dblarr(5), $
                   dblarr(5), dblarr(5), dblarr(5), dblarr(5) )

out = replicate(oe, n_elements(d))
struct_assign, d, out

for i=0,n_elements(d)-1 do begin
    print, i
    out[i].TOTAL_FLUX[0] = totalsersicflux(out[i].BULGEFIT,cutoff=0)
    out[i].TOTAL_FLUX[1] = totalsersicflux(out[i].DISKFIT,cutoff=0)
    out[i].TOTAL_FLUX[2] = totalsersicflux(out[i].DVCFIT,cutoff=0)
    out[i].TOTAL_FLUX[3] = totalsersicflux(out[i].EXPFIT,cutoff=0)
    out[i].TOTAL_FLUX[4] = totalsersicflux(out[i].SERSICFIT,cutoff=0)

    out[i].TOTAL_FLUX_ERR[0] = totalsersicflux_error(out[i].BULGEFIT,out[i].COVAR_BULGE,cutoff=0)
    out[i].TOTAL_FLUX_ERR[1] = totalsersicflux_error(out[i].DISKFIT,out[i].COVAR_DISK,cutoff=0)
    out[i].TOTAL_FLUX_ERR[2] = totalsersicflux_error(out[i].DVCFIT,out[i].COVAR_DVC,cutoff=0)
    out[i].TOTAL_FLUX_ERR[3] = totalsersicflux_error(out[i].EXPFIT,out[i].COVAR_EXP,cutoff=0)
    out[i].TOTAL_FLUX_ERR[4] = totalsersicflux_error(out[i].SERSICFIT,out[i].COVAR_SERSIC,cutoff=0)

    out[i].BTT[0] = bulgetotot(out[i].BULGEFIT,cutoff=0)
    out[i].BTT[1] = bulgetotot(out[i].DISKFIT,cutoff=0)
    out[i].BTT[0] = bulgetotot_error(out[i].BULGEFIT, $
                                     out[i].COVAR_BULGE)
    out[i].BTT[1] = bulgetotot_error(out[i].DISKFIT,$
                                     out[i].COVAR_DISK)


    out[i].BULGE_FLUX = out[i].TOTAL_FLUX*out[i].BTT
    out[i].DISK_FLUX = out[i].TOTAL_FLUX*(1-out[i].BTT)
    
    out[i].BULGE_FLUX_ERR[0] = totalsersicflux_error(out[i].BULGEFIT,out[i].COVAR_BULGE,cutoff=0,/bulge)
    out[i].BULGE_FLUX_ERR[1] = totalsersicflux_error(out[i].DISKFIT,out[i].COVAR_DISK,cutoff=0,/bulge)
    out[i].DISK_FLUX_ERR[0] = totalsersicflux_error(out[i].BULGEFIT,out[i].COVAR_BULGE,cutoff=0,/disk)
    out[i].DISK_FLUX_ERR[1] = totalsersicflux_error(out[i].DISKFIT,out[i].COVAR_DISK,cutoff=0,/disk)
    
endfor

mwrfits, out, file_out, /create

end
