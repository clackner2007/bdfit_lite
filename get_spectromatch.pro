;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get matching object from VAGC spectro list
;
; all objects should be in that list
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO get_spectromatch, in_file, spec_file, out_file

spectro=mrdfits(spec_file,1,columns=['SDSS_SPECTRO_TAG','PLUG_RA',$
                                     'PLUG_DEC', 'PLATE', 'MJD',$
                                     'FIBERID', 'Z', 'Z_ERR'])

mydata = mrdfits(in_file,1,columns=['RA','DEC', 'Z'])
num = n_elements(mydata)

;match the lists, use 3" as max length
spherematch, mydata.RA, mydata.DEC, spectro.PLUG_RA, $
  spectro.PLUG_DEC, 4*8.4e-4, mydataMatch, specMatch

;print, mydataMatch
;print, specMatch
;print, specMatch[sort(mydataMatch)]
if n_elements(mydataMatch) ne num then print, $
  "didn't match all the objects", num,' matched ',n_elements(mydataMatch)

ss=sort(mydataMatch)
for n=0L, n_elements(mydata)-1 do begin
    xx=where(mydataMatch[ss] eq n, nx)
    if nx eq 0 then $
      print, n, where(mydataMatch eq n)
endfor

mwrfits, spectro[specMatch[ss]], out_file, /create


END
