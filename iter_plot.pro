;;;;;;;;;;;;;;;;;;;;;;;;;
; iter_plot.pro
; plots oned profile of npfit iteration
;;;;;;;;;;;;;;;;;;;;;;;;;

PRO iter_plot, myfunct, p, iter, fnorm, FUNCTARGS=functargs, $
               parinfo=parinfo, quiet=quiet,_EXTRA=_EXTRA
print, 'iteration: ',iter,' chi-squared: ',fnorm
print, 'params: ',p
print, '+++++++++++++++++++++++++++++++++'

if( n_elements(p) eq 14 ) then begin
    makegrid, 2*p[12], 2*p[13], xx,yy
    psf= _EXTRA.psf
    oned_image, flux, rad, iv, pixelfluxpsf(xx,yy,p,_EXTRA={cutoff:1,psfImage:psf}), $
      p[11],p[9]
    oplot, alog10(rad+1.0), flux
endif
END
