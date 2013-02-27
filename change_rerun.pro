;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; use spherematch to change between reruns
;
;  returns new rerun, objid list
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO change_rerun, run, camcol, field, ra, dec, r_mag, rerun=rerun,$
                  objid=objid, properties=properties, tags=tags


BRIGHT=2L^01
BLENDED=2L^03
NODEBLEND=2L^06
PRIMARY = sdss_flagval('RESOLVE_STATUS','RUN_PRIMARY')
DEBLEND_AS_PSF = 2L^25
STAR=6
CR=1
GALAXY=3

thresh = 2.0/3600.

if keyword_set(tags) then begin
  junk=where(tags eq 'NCHILD', nc)
  f=where(tags eq 'OBJC_FLAGS', nf)
  r=where(tags eq 'RESOLVE_STATUS', nr)
  gs=where(tags eq 'OBJC_TYPE', ng)
  bar=where(tags eq 'MODELFLUX', nm)
  if nc eq 0 then tags=[[tags],'NCHILD']
  if nf eq 0 then tags=[[tags],'OBJC_FLAGS']
  if nr eq 0 then tags=[[tags],'RESOLVE_STATUS']
  if ng eq 0 then tags=[[tags], 'OBJC_TYPE']
  if nm eq 0 then tags=[[tags], 'MODELFLUX']
  stags=['RA','DEC','RERUN','ID',tags]
endif else $
  stags=['RA','DEC','RERUN','ID', 'NCHILD', 'OBJC_FLAGS', $
         'RESOLVE_STATUS', 'OBJC_TYPE','MODELFLUX']


if n_elements(run) eq 1 then begin
    new = sdss_readobj(run, camcol, field, rerun=rerun, $
                       select_tags=stags)
    spherematch, ra, dec, new.RA, new.DEC, thresh, m1, m2
    objid=new[m2].ID
    properties = new[m2]
    
endif else begin

    objid = intarr(n_elements(run)) - 1
    runs = SetIntersection(run, run)
    prop_e = -1
    id = 0
    while not is_struct(prop_e) and id lt n_elements(run) do begin
        prop_e = sdss_readobj(run[id], camcol[id], field[id], $
                              rerun=rerun, select_tags=stags)
        id += 1
    endwhile
    if not is_struct(prop_e) then begin
        print, "no objects found"
        return
    endif
    properties = replicate(prop_e[0], n_elements(objid))
    properties.NCHILD = -1
    
    for i=0L, n_elements(runs)-1 do begin
        for col = 1, 6 do begin
            current=where(run eq runs[i] and camcol eq col, ncr)
            if ncr lt 1 then continue

          fields = SetIntersection(field[current], field[current])

        ;    for f=0, n_elements(fields)-1 do begin

          ;      current = where(run eq runs[i] and camcol eq col and $
           ;                     field eq fields[f])
            print, "run ",runs[i]," camcol ",col
            new=sdss_readobj(runs[i], col, $
                             fields, $
                             rerun=rerun, $
                             select_tags=stags, $
                             /silent )

            if not is_struct(new) then continue
            ;don't count bright (r) objects or objects with NCHILD !=
            ;0
            ;ignore objects which are BLENDED and don't have NODEBLEND
            keep = where((new.NCHILD eq 0) and $
                         ((new.OBJC_FLAGS and BRIGHT) eq 0) and $
                         (((new.OBJC_FLAGS and BLENDED) eq 0) or $
                          ((new.OBJC_FLAGS and NODEBLEND) ne 0)) and $ 
                         ((new.RESOLVE_STATUS and PRIMARY) ne 0) and $
                         ((new.OBJC_TYPE eq GALAXY)), nk) ;and $
                         ;((new.OBJC_FLAGS and DEBLEND_AS_PSF) eq 0),$
                         ;nk)
            if nk le 0 then continue
            new = new[keep]
            
            spherematch, ra[current], dec[current], new.RA, new.DEC, $
              thresh, m1, m2, maxmatch=0
            if m1[0] ne -1 then begin
                for mm=0L, n_elements(current) do begin
                    matches = where(m1 eq mm,nmm)
                    if nmm eq 0 then continue
                    if n_elements(matches) eq 1 then begin
                        objid[current[m1[matches]]] = new[m2[matches]].ID
                        properties[current[m1[matches]]] = $
                          new[m2[matches]]
                    endif else begin
                        magdiff = abs(r_mag[current[m1[matches]]]$
                                      -new[m2[matches]].MODELFLUX[2])
                        print, magdiff
                        print, sdss_galid(ra[current[m1[matches[0]]]],$
                                          dec[current[m1[matches[0]]]])
                        closest = $
                          where(abs(magdiff - min(magdiff)) lt 1.e-9)
                        objid[current[m1[matches[closest]]]] = new[m2[matches[closest]]].ID
                        properties[current[m1[matches[closest]]]] = $
                          new[m2[matches[closest]]]
                    endelse    
                endfor
                
            endif
        endfor
        
    endfor
endelse




END
