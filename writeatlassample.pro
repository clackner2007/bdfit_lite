;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
; Jan 11 2010
; 
;  write the atlas sample to a set of files
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO writeatlassample, output


chunk = 20000L
number = 0L*20000L
max=2506754L

while (number le max) do begin
    print, number,' of 2.5e6'
    last = min([number+chunk, 2506754L])
    list = atlassample_extra([number, last], [0.003,0.05], rerun=301)
    if not is_struct(list) then begin
        number += chunk
        continue
    endif
    size = n_elements(list)
    name = output + string(number/chunk, size, $
                           format='(i04,"-",i05,".fits")')
    number += chunk
    mwrfits, list, name
endwhile

END
