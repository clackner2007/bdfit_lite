;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get matching object from VAGC list and get Two mass kband data
;
; and add them to in_file data, and dump to out_file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO get_twomass_match, in_file, vagc_file, twomass_file, out_file

vagc=mrdfits(vagc_file,1,columns=['RA','DEC','TWOMASS_TAG','TWOMASS_POSITION'])

mydata = mrdfits(in_file,1);,columns=['RA','DEC','RUN','RERUN',$
                           ;         'CAMCOL','FIELD','ID','Z',$
                           ;        "PETROFLUX",$
                           ;         "PETROFLUX_IVAR",$
                           ;         "MODELFLUX",$
                           ;         "MODELFLUX_IVAR",$
                           ;         "PSFFLUX",$
                           ;         "PSFFLUX_IVAR",$
                           ;         "EXTINCTION"])
num = n_elements(mydata)

;match the lists, use 3" as max length
spherematch, mydata.RA, mydata.DEC, vagc.RA, vagc.DEC, 3*8.4e-4, mydataMatch, vagcMatch


if n_elements(mydataMatch) ne num then $
  print, "didn't match all the objects, from ", num,' matched ',n_elements(mydataMatch)

twomass = mrdfits(twomass_file,1,columns=["RA", $
                                          "DECL", $
                                          "J_M_EXT",$
                                          "J_MSIG_EXT",$
                                          "J_FLG_EXT",$
                                          "H_M_EXT",$
                                          "H_MSIG_EXT",$
                                          "H_FLG_EXT",$
                                          "K_M_EXT",$
                                          "K_MSIG_EXT",$
                                          "K_FLG_EXT"])



output_e =create_struct(mydata[0], ["TMASS_K_M", $
                                    "K_KCORRECT",$
                                    "TMASS_J_M",$
                                    "J_KCORRECT",$
                                    "TMASS_H_M",$
                                    "H_KCORRECT"],$
                        0.,0.,0.,0.,0.,0.)
                        
output = replicate(output_e,n_elements(mydataMatch))
struct_assign, mydata[mydataMatch], output



;output.RUN = mydata[mydataMatch].RUN
;output.RERUN = mydata[mydataMatch].RERUN
;output.CAMCOL = mydata[mydataMatch].CAMCOL
;output.FIELD = mydata[mydataMatch].FIELD
;output.ID = mydata[mydataMatch].ID
;output.Z = mydata[mydataMatch].Z
;output.FITOUT_ID = mydataMatch

hasTM = where(vagc[vagcMatch].TWOMASS_POSITION gt 0)
;positions = (vagc[vagcMatch].TWOMASS_POSITION)[hasTM]
positions = (vagcMatch[hasTM])


;if missing sdss, go get it
if (tag_indx(mydata, 'PETROFLUX') eq -1) then begin
    my_sdss_data = $
      mrdfits('/peyton/scr/depot1/clackner/vagc2/object_sdss_imaging.fits',$
              1,rows=positions)

    Kcorrect = Transpose((twomass_kcorrect(output[hasTM].Z,$
                                           twomass=twomass[positions],$
                                           calibobj=my_sdss_data)))
   
endif else begin
    Kcorrect = Transpose((twomass_kcorrect(output[hasTM].Z,$
                                       twomass=twomass[positions],$
                                       calibobj=(mydata[mydataMatch])[hasTM])))
endelse

;report AB magnitudes
output[hasTM].TMASS_K_M = twomass[positions].K_M_EXT+1.84
output[hasTM].TMASS_J_M = twomass[positions].J_M_EXT+0.894
output[hasTM].TMASS_H_M = twomass[positions].H_M_EXT+1.374

sdss_bands=5
output[hasTM].K_KCORRECT = Kcorrect[*,7-sdss_bands]
output[hasTM].J_KCORRECT = Kcorrect[*,5-sdss_bands]
output[hasTM].H_KCORRECT = Kcorrect[*,6-sdss_bands]

ss = sort(mydataMatch)
output = output[ss]


mwrfits, output, out_file, /create



END
