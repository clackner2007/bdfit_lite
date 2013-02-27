;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; extracts a sample from the total atlassample2
; which is constrained by in absolute magnitude and apparent magnitude
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO get_mag_sample, inputfolder, inputfile, absLow, absHigh, $
                    apprLow, apprHigh, nfiles=nfiles, band=band

files = unixfind(inputfolder, inputfile )
hubble = 70./3.0e5/1.0e5

if not keyword_set(band) then band=2
if not keyword_set(nfiles) then nfiles=1e5

print, "number of files ",n_elements(files)
print, files
for i=0, min([n_elements(files)-1,nfiles]) do begin
    temp = mrdfits(files[i], 1, /silent)
    apprmag = -2.5*alog10(temp.MODELFLUX[band]*1.0e-9)
    distmod = 5.0*alog10(lumdist(temp.Z, /silent)*1.0e5)
    absMag = apprmag - temp.EXTINCTION[band] - distmod
    
    meets_crit = where(absMag ge absLow and $
                       absMag lt absHigh and $
                       apprmag ge apprLow and $
                       apprmag lt apprHigh)
    if(meets_crit[0] eq -1 ) then continue
    if n_elements(data) eq 0 then data = temp[meets_crit] $
    else data=[[data],temp[meets_crit]]
     
endfor

print, 'number found: ',n_elements(data)
outputname = inputfolder+string(apprLow, apprHigh, abs(absLow), abs(absHigh), $
                                format='("m",i02,"_",i02,"M-",i02,"_",i02,".fits")')

mwrfits, data, outputname, /create

END
