;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get matching object from VAGC photo list
;
; all objects should be in that list
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO get_photomatch, in_file, photo_file, param_file, out_file

photo=mrdfits(photo_file,1)
mydata=mrdfits(in_file,1)

num = n_elements(mydata)

;keep
keep=where((photo.NCHILD eq 0))


;match the lists, use 3" as max length
spherematch, mydata.RA, mydata.DEC, photo[keep].RA, $
  photo[keep].DEC, 2.4e-4, mydataMatch, photoMatch

;print, mydataMatch
;print, specMatch
;print, specMatch[sort(mydataMatch)]
if n_elements(mydataMatch) ne num then begin
    print, $
      "didn't match all the objects from", num,' matched ',n_elements(mydataMatch)
    return
endif
ss=sort(mydataMatch)
for n=0L, n_elements(mydata)-1 do begin
    xx=where(mydataMatch[ss] eq n, nx)
    if nx eq 0 then $
      print, n, where(mydataMatch eq n)
endfor


objs=(photo[keep])[photoMatch[ss]]

e=create_struct(['GROUP_ID', "RA", "DEC", 'Z', "FIBER_COL",$
                 "RUN", "CAMCOL", "FIELD", "ID", "RERUN", $
                 "ABSMAG", "KCORRECT_0"], 0L, 0.0D, 0.0D, $
                0.0, 0, 0, 0, 0, 0, 0, dblarr(5), dblarr(5))

output = replicate(e, n_elements(objs))
output.GROUP_ID = mydata.GROUP_ID
output.RA = mydata.RA
output.DEC = mydata.DEC
output.Z = mydata.Z
output.FIBER_COL = mydata.FIBER_COL
output.RUN = objs.RUN
output.CAMCOL = objs.CAMCOL
output.FIELD = objs.FIELD
output.ID = objs.ID
output.RERUN = objs.RERUN

distMod = 5.0*alog10(lumdist(output.Z, /silent)*1.e5)
darray=dblarr(5,num)
b=[1.4e-10,0.9e-10,1.2e-10,1.8e-10,7.4e-10]
barray= dblarr(5,num)+1.0
for i=0, 4 do begin
    darray[i,*]=distMod
    barray[i,*] *= b[i]
endfor

output.KCORRECT_0 = sdss_kcorrect(output.Z, calibobj=objs)
output.ABSMAG = -2.5*alog10(objs.MODELFLUX*1.e-9) - $
  output.KCORRECT_0 - objs.EXTINCTION - darray
;(asinh((objs.MODELFLUX*1.e-9/(2.0*barray)) + alog(barray)) $
                ; *(-2.5/alog(10.0)) - output.KCORRECT_0 $
                 ;- output.EXTINCTION - darray)


mwrfits, output, out_file, /create

END
