;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
; Mar 3 2010
; plotatlasimages.pro
;
;
; plots the atlas images, as well as their B+D decompositions
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO plotatlasimages, inputfile, filter=filter, residuals=residuals

;default to r-band
if n_elements(filter) eq 0 then filter=2

data = mrdfits( inputfile, 1, /silent )

!P.MULTI = [0,3,n_elements(data)]
!P.CHARSIZE=1.8
;!Y.OMARGIN = [1,3]
;!X.OMARGIN = [2,2]
;    set_plot, "x"
set_plot, "ps"
bands=['u','g','r','i','z']
last = fix((strsplit(inputfile,'.',/extract))[1])
output = string(last-9,last,bands[filter],format='("atlassampleDR7/",i04,".",i04,"-",a,".eps")')
if( keyword_set(residuals) ) then begin ;plot residuals instead of
                                        ;images
    output=string(last-9,last,bands[filter],$
                  format='("atlassampleDR7/",i04,".",i04,"-",a,"-resid.eps")')
endif
print, 'output to: ',output
device, /encapsulated, filename=output, /inches
device, bits_per_pixel=8

for i=0L, n_elements(data) - 1 do begin
    atlas = (sdss_atlas_image( data[i].RUN, $
                              data[i].CAMCOL, $
                              data[i].FIELD, $
                              data[i].ID, $
                              rerun=data[i].RERUN, $
                              filter=filter ))*1.0D
    psffile = sdss_name( 'psField', data[i].RUN, $
                         data[i].CAMCOL, data[i].FIELD, $
                         rerun=data[i].RERUN )
    locpsf = sdss_psf_recon(mrdfits(psffile, filter+1, /silent), $
                            data[i].ROWC+0.5, $
                            data[i].COLC+0.5, normalize=1.0 )
    
    lenx = (size(atlas, /dimensions))[0]
    leny = (size(atlas, /dimensions))[1]
    if( i eq 0 ) then device, xsize=20.0, ysize=60.

    makegrid, (size(atlas, /dimensions))[0], $
      (size(atlas, /dimensions))[1], $
      x, y

;     n4bulge = [data[i].BULGEFIT[6:10, filter], $
;                data[i].BULGEFIT[12:13, filter], $
;                data[i].BULGEFIT[11, filter] ]
;     n4disk = [data[i].BULGEFIT[0:4, filter], $
;               data[i].BULGEFIT[12:13, filter], $
;               data[i].BULGEFIT[5, filter] ]
;     n1bulge = [data[i].DISKFIT[6:10, filter], $
;                data[i].DISKFIT[12:13, filter], $
;                data[i].DISKFIT[11, filter] ]
;     n1disk = [data[i].DISKFIT[0:4, filter], $
;               data[i].DISKFIT[12:13, filter], $
;               data[i].DISKFIT[5, filter] ]

    if( data[i].BULGEFIT[1,filter] lt 1.0e-9  or $
        data[i].BULGEFIT[7, filter] lt 1.0e-9 ) then begin
        bulgedisk = dblarr((size(atlas, /dimensions))[0], $
                           (size(atlas, /dimensions))[1] )+1000.0
        bulgeToTotn4 = -1.0
    endif else begin
        bulgedisk = pixelfluxpsf( x, y, data[i].BULGEFIT[*,filter], $
                                  _EXTRA={psfImage:locpsf} )+1000.0
;         bulgeToTotn4 = totalsersicflux(n4bulge)/$
;           (totalsersicflux(n4bulge) + totalsersicflux(n4disk))
        bulgeToTotn4=bulgetotot(data[i].BULGEFIT[*,filter])
    endelse 
    if( data[i].DISKFIT[0,filter] lt 1.0e-9 or $
        data[i].DISKFIT[7, filter] lt 1.0e-9 ) then begin ;or $
  ;      data[i].DISKFIT[6,filter] gt 1.0e-3 ) then begin
        diskdisk = dblarr((size(atlas, /dimensions))[0], $
                          (size(atlas, /dimensions))[1] )+1000.0
        bulgeToTotn1 = -1.0
    endif else begin
        diskdisk = pixelfluxpsf( x, y, data[i].DISKFIT[*,filter], $
                                 _EXTRA={psfImage:locpsf} )+1000.0
;         bulgeToTotn1 = totalsersicflux(n1bulge)/$
;           (totalsersicflux(n1bulge) + totalsersicflux(n1disk))
        bulgeToTotn1 = bulgetotot(data[i].DISKFIT[*,filter])
    endelse

    chi = '!4'+String("166B)+'!X'
    
    imagename = string(data[i].RUN, data[i].CAMCOL, data[i].FIELD, $
                       data[i].ID,format='(i04,".",i0,".",i04,".",i04," ")')
    name2=string(chi, data[i].CHISQ_BULGE[filter],bulgeToTotn4,$
                 format='("dVc ",a0,"^2=",f0.2,",b/T=",f0.2," ")')
    name3=string(chi,data[i].CHISQ_DISK[filter],bulgeToTotn1,$
                 format='("exp ",a0,"^2=",f0.2,",b/T=",f0.2)')

    
    !P.TITLE=imagename
    tvscale, atlas, /keep_aspect_ratio, multimargin=[1,0,2,0]
;,position=[0.0,0.0,0.32,0.9]
    !P.TITLE=name2
    if( keyword_set(residuals) ) then $
      tvscale, bulgedisk-atlas, /keep_aspect_ratio, $
      multimargin=[1,0,2,0] $
    else $
      tvscale, bulgedisk, /keep_aspect_ratio, multimargin=[1,0,2,0]
;, position=[0.34,0.0,0.66,0.9]
    !P.TITLE=name3
    if(keyword_set(residuals) ) then $
      tvscale, diskdisk-atlas, /keep_aspect_ratio, $
      multimargin=[1,0,2,0] $
    else $
      tvscale, diskdisk, /keep_aspect_ratio, multimargin=[1,0,2,0]
;,position=[0.68,0.0,1.0,0.9]
;     xyouts, 0.05, 0.92, alignment=0.0, charsize=1., /normal,imagename
;     xyouts, 0.35, 0.92, alignment=0.0, charsize=1.,/normal,name2
;     xyouts,0.7,0.92,alignment=0.0,charsize=1.,/normal,name3

;    set_plot, "x"
endfor

device, /close
!P.MULTI=0
!Y.OMARGIN = [0,0]
!P.TITLE=''
!P.CHARSIZE=1.0
END
