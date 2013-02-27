;____________________________________
;
; Claire Lackner
;  Nov, 5, 2009
;
;  readoutput.pro
;
;  IDL
; reads the output tables from the fitter
;   
;____________________________________


;uses gals which is list of galaxies from sdss lowz catalog
;takes 4 arguments, the ranges in which to look for sources in the
;lowz catalog, both given in degrees
PRO readfakes, folder, data

;usage statement
if( n_params() ne 2 ) then message, $
  'Usage: readoutput, folder_of_output_files'

filenames = unixfind(folder, "fake*.fits")

temp = mrdfits( filenames[0], 1 )
data = temp
;sersic_n = temp.SERSIC_N
;blanton_n = temp.BLANTON_N

for i=1L, n_elements( filenames )-1 do begin
    temp = mrdfits( filenames[i], 1)
    data = [[data], temp]
;    sersic_n = [[sersic_n], [temp.SERSIC_N]]
;    blanton_n = [[blanton_n], [temp.BLANTON_N]]
endfor

sersic = where( data.TYPE eq 'sersic    '  and data.SERSIC_N[0] > 1.0e-2)
bulges = where( data.TYPE eq 'bulgefdisk' and data.SERSIC_N[0] > 1.0e-2)
pseudo = where( data.TYPE eq 'fdiskfdisk' and data.SERSIC_N[0] > 1.0e-2)

tvlct, 0,0,0,0
tvlct, 255, 0, 0, 1 ;red
tvlct, 0,255,0,2 ;green
tvlct,0,0,255,3 ;blue
tvlct,255,150,0,4 ;orange
tvlct,255,255,0,5;yellow
tvlct,255,0,255,6;purple
tvlct,0,255,255,7;light b
tvlct,255,170,170,8;pink
tvlct,0,90,30,9;dark green
tvlct,130,90,40,10;brown

set_plot, "ps"
device, /encapsulated, /color, filename='out.eps'

plot, data[sersic].BULGEPARAMETERS[2], data[sersic].SERSIC_N[0], $
  psym=1, color=0, xrange=[0.1,5.], yrange=[0.1,17.]
oplot, data[sersic].BULGEPARAMETERS[2], data[sersic].SERSIC_N[1],$
psym=2, color=1
oplot, data[sersic].BULGEPARAMETERS[2], data[sersic].SERSIC_N[2],$
psym=2, color=2
oplot, data[sersic].BULGEPARAMETERS[2], data[sersic].SERSIC_N[3],$
psym=2, color=3
oplot, data[sersic].BULGEPARAMETERS[2], data[sersic].SERSIC_N[4],$
psym=2, color=4
oplot, data[sersic].BULGEPARAMETERS[2], data[sersic].SERSIC_N[5],$
psym=2, color=5
oplot, data[sersic].BULGEPARAMETERS[2], data[sersic].SERSIC_N[6],$
psym=2, color=6
oplot, data[sersic].BULGEPARAMETERS[2], data[sersic].SERSIC_N[7],$
psym=2, color=7
oplot, data[sersic].BULGEPARAMETERS[2], data[sersic].SERSIC_N[8],$
psym=2, color=9

device, /close
set_plot, "x"

; print, size( blanton_n, /dimensions)
; x = [0,15]
; y=x
; set_plot, "ps"
; device, /encapsulated, filename="sersic_n_iband.eps"

; plot, sersic_n[3, *], blanton_n[3, *], psym=1, title="sersic index, i band",  $
;   xtitle="two-d fitter, constant weights", ytitle="blanton, 2003"
; oplot, x,y
; device, /close

; set_plot, "x"

; good = where( data.SERSIC_N[2] GT 0 )

; set_plot, "ps"
; device, /encapsulated, /color, filename="sersic_n_ranges.eps"

;plot a histogram of sersic indices
; plot, indgen(20)*9.0/20.0 + 0.1, $
;   histogram(data[good].sersic_n[2], min=0.1, max=8.0, nbins=20 ), $
;   xtitle='Sersic Index', ytitle='Number of galaxies'
; oplot, indgen(20)*9.0/20.0 + 0.1, $
;   histogram(data[good].fs_sersic_n[1,2], min=0.1, max=8.0, nbins=20),$
;   linestyle=1,color=1
; oplot, indgen(20)*9.0/20.0 + 0.1, $
;   histogram(data[good].fs_sersic_n[2,2], min=0.1, max=8.0, nbins=20),$
;   linestyle=2,color=2
; oplot, indgen(20)*9.0/20.0 + 0.1, $
;   histogram(data[good].fs_sersic_n[3,2], min=0.1, max=8.0, nbins=20),$
;   linestyle=3, color=3

;print, max(sersic_n[3,*]), min(sersic_n[3,*])

; device, /close
; set_plot, "x"

END

PRO plottype, type, key1, key2, data, bulge, disk

tag1 = (where( tag_names(data) eq key1 ))[0]
tag2 = (where( tag_names(data) eq key2 ))[0]

want = where(data.TYPE eq type and data.NORM[0] gt 1.0e-4 )

tvlct, 0,0,0,0
tvlct, 255, 0, 0, 1 ;red
tvlct, 0,255,0,2 ;green
tvlct,0,0,255,3 ;blue
tvlct,255,150,0,4 ;orange
tvlct,255,255,0,5;yellow
tvlct,255,0,255,6;purple
tvlct,0,255,255,7;light b
tvlct,255,170,170,8;pink
tvlct,0,90,30,9;dark green
tvlct,130,90,40,10;brown

if type eq 'sersic    ' then type='sersic'
set_plot, "ps"
file=key1+'vs'+key2+type+'.eps'
device, /encapsulated, /color, filename=file

if bulge[0] ne -1 then begin
    b = data[want].BULGEPARAMETERS[bulge]
endif
if disk[0] ne -1 then begin
    d=data[want].DISKPARAMETERS[disk]
endif

data1 = data[want].(tag1)
data2 = data[want].(tag2)

xmax = max(data1) +0.5
xmin = min(data1) - 0.5
ymax = max(data2) + 0.5
ymin = min(data2) - 0.5

plot, data1[0,*], data2[0,*], color=0, psym=1, xrange=[xmin,xmax], $
  yrange=[ymin,ymax], xtitle=key1, ytitle=key2, title=file

for i=1, 8 do begin
    oplot, data1[i,*], data2[i,*], color=i+1, psym=1
endfor

if bulge[0] ne -1 then begin
    oplot, b[0,*], b[1,*], color=1, linestyle=0
endif
if disk[0] ne -1 then begin
    oplot, d[0,*], d[1,*], color=1, linestyle=1
endif

device, /close
set_plot, "x"

END
