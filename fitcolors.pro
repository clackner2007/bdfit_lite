;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
; Apr 6 2010
;
; fitcolors.pro
;
;
; given a fit in one band, fit in separate bands by simply scaling the
; bulge and the disk fluxes. Return the new fit normalizations, new B/T in
; each band
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION fitcolors, object, basisparams, bands

newparams = dblarr(n_elements(basisparams), n_elements(bands))

for i=0, n_elements(bands)-1 do begin
    atlas = (sdss_atlas_image( object.RUN, object.CAMCOL, $
                               object.FIELD, object.ID, $
                               rerun=object.RERUN, $
                               filter = bands[i] ))*1.0D
    psffile = sdss_name('psField', object.RUN, object.CAMCOL, $
                        object.FIELD, rerun=object.RERUN)
    psf = sdss_psf_recon(mrdfits(psffile,bands[i]+1, /silent), $
                            object.ROWC+0.5, object.COLC+0.5, $
                            normalize=1.0D)
    ivar = get_atlas_ivar( atlas, object.RUN, object.CAMCOL, $
                           object.FIELD, object.ID, $
                           rerun=object.RERUN, filter=bands[i] )
    im_unbias = (atlas-1000.)*1.0D
    im_size = size(atlas, /dimensions)

    ;use the linear method only
;    newparams[*,i] = scalefitMP( im_unbias, psf, ivar, basisparams )

    newparams[*,i] = scalefit( im_unbias, psf, ivar, basisparams )

endfor
return, newparams
END
