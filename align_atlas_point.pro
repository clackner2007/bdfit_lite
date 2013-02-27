;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; align_image.pro
;
; Claire Lackner
; Aug 13 2010
; code to get atlas image alignments from fpC files
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO align_atlas_point, x, y, x_out, y_out, filter_in, filter_out, $
                       run, camcol, field, id, rerun=rerun,$
                       fpc_file=fpc_file, fpcheader=fpcheader

;at=sdss_atlas_construct()
at = sdss_atlas(sdss_name('fpAtlas', run, camcol, field, rerun=rerun), id)

cmin = at.MASTER_CMIN+at.DCOL[[filter_in,filter_out]]
rmin = at.MASTER_RMIN+at.DROW[[filter_in,filter_out]]
sdss_atlas_destroy, at

if keyword_set(fpcheader) then begin
    fpC_in = *(fpcheader[filter_in])
    fpC_out = *(fpcheader[filter_out])
endif else begin
    if keyword_set(fpc_file) then begin
        data=mrdfits(fpc_file,filter_in,hdr,/silent)
        fpC_in = hdr
        data=mrdfits(fpc_File,filter_out,hdr,/silent)
        fpC_out = hdr
    endif else begin
        ims= sdss_readfpc(run,camcol,field,rerun=rerun,filter=[filter_in,filter_out], $
                          phdr=p_hdr)
        
        fpC_in = *(p_hdr[0])
        fpC_out = *(p_hdr[1])
        heap_free, p_hdr[0]
        heap_free, p_hdr[1]
    endelse
endelse


;difference from fixed RA/DEC point in input band to point of interest
x_atlas_in = x+cmin[0]
y_atlas_in = y+rmin[0]

xyxy, fpC_in, fpC_out, x_atlas_in, y_atlas_in, x_atlas_out, y_atlas_out

x_out=x_atlas_out - cmin[1]
y_out=y_atlas_out - rmin[1]


END
