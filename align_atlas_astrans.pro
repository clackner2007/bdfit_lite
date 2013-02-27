;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; align_atlas_astrans.pro
;
; Claire Lackner
; Aug 13 2010
; aligns atlas images using astrans structures
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO align_atlas_astrans, x, y, x_out, y_out, filter_in, filter_out, $
                       run, camcol, field, id, rerun=rerun

;at=sdss_atlas_construct()
at = sdss_atlas(sdss_name('fpAtlas', run, camcol, field, rerun=rerun), id)

cmin = at.MASTER_CMIN+at.DCOL[[filter_in,filter_out]]
rmin = at.MASTER_RMIN+at.DROW[[filter_in,filter_out]]
sdss_atlas_destroy, at


;difference from fixed RA/DEC point in input band to point of interest
x_atlas_in = x+cmin[0] -0.5 ;sdss not zero-indexed
y_atlas_in = y+rmin[0] -0.5


;get the exact astrometry for the asTrans file
astrans_in = sdss_exactastrom(run,camcol,field,$
                              rerun=rerun,filter=filter_in)
astrans_out = sdss_exactastrom(run,camcol,field, $
                              rerun=rerun,filter=filter_out)

astrans_xy2eq, astrans_in, x_atlas_in, y_atlas_in, $
  ra=ra, dec=dec, mu=mu, nu=nu

astrans_eq2xy, astrans_out, ra, dec, xpix=x_atlas_out, $
  ypix=y_atlas_out, mu=mu, nu=nu

;move to atlas image coordinates and un-zero index coordinates
x_out=x_atlas_out - cmin[1]+0.5
y_out=y_atlas_out - rmin[1]+0.5


END
