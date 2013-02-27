;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
;
; add_kcorrect.pro
;
;  adds the kcorrection to the processed data
;  kcorrects each band to z=0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO add_kcorrect, procfile, rawfile

data = mrdfits(procfile, 1)
calibobjs = mrdfits(rawfile,1)

new_entry = create_struct(data[0], ['KCORRECT_0'], dblarr(5))

table = replicate(new_entry, n_elements(data))
struct_assign, data, table

;;; calculate the k correction using KCORRECT

table.KCORRECT_0 = sdss_kcorrect(data.Z, calibobj=calibobjs)

mwrfits, table, procfile, /create

END
