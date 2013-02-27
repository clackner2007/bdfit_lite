;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
; Apr 13 2010
;
;  sdss_galid.pro
;
;  ARGUMENTS: RA (J2000 in degrees)
;             DEC (J200 in degrees)
;
;  OPTIONAL KEYWORDS: SHORT
;     set this to return the shortend name: SDSSJhhmmss[+-]dd
;
;  RETURN: IAU name of galaxy in form:
;    SDSSJhhmmss.ss[+-]ddmmss.s
;   where the first number is the RA, and the second is the DEC
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



function sdss_galid, ra, dec, short=short

;the right ascension is always positve
hrRA = floor(RA/360.*24.)
minRA = floor((hrRA - (RA/360.*24.))*60.*(-1.))
secRA = truncd(((RA/360.*24. - hrRA)*60. - minRA)*60., (-2))

;need to worry about signs for dec
degDEC = floor( abs(DEC) )*sgn( DEC )
minDEC = floor(abs(DEC - degDEC)*60.)
secDEC = truncd((abs(DEC-degDEC)*60. - minDEC)*60., (-1))

posneg = string( sgn(DEC) gt 0? '+': '-' )

if( keyword_set(short) ) then $
  name = string(hrRA, minRA, secRA, posneg, abs(degDEC), $
              format='("SDSSJ",i02,i02,i02,a,i02)') $
  else $
  name = string(hrRA, minRA, secRA, posneg, abs(degDEC), $
                minDEC, secDEC, $
                format='("SDSSJ",i02,i02,f05.2,a,i02,i02,f04.1)' )

return, name

end
