;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
; Feb 28 2010
;
; get_atlas_ivar.pro
;
;
; returns the inverse variance image of an atals image, the formula
; used to calculate the varaiance is:
; var = darkVarainace + (objCounts + skyCounts)/gain
; the darkVariance and gain are in the header of the fpObjc file, the
; sky counts are in the fpObjc file, stored for each object
;
;Currently, this function does not account for pixels that have been
;interpolated over (which have different noise properties), but that
;effect should be small, especially for bright galaxies
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;get the sky error by calculating the sky in the neighboring 128x128
;squares in the frame and taking the variance of the neighbors, this
;is a rough estimate, but eliminates the ridiculously large sky_err
;values 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION get_128sky_var, run, camcol, field, rowc, colc, $
                         rerun=rerun, filter=filter

skyvals=dblarr(8)
nvals=0

;get the sky from the second HDU of the fpBIN file, value is shifted
;by 10 bits
sky = mrdfits(sdss_name('fpBIN', run, camcol, field, $
                        rerun=rerun,filter=filter)+".gz",$
              1, /silent)/1024.

xi = floor(colc/128.)
yi = floor(rowc/128.)
;print, rowc, colc, xi, yi

maxX = (size(sky, /dimension))[0]-1
maxY = (size(sky, /dimension))[1]-1

if xi ne 0 then begin
    skyvals[nvals] = sky[xi-1, yi]
    nvals += 1
    if yi ne 0 then begin
        skyvals[nvals] = sky[xi-1,yi-1]
        nvals+=1
    endif
    if yi ne maxY then begin
        skyvals[nvals] = sky[xi-1, yi+1]
        nvals += 1
    endif
endif
if xi ne maxX then begin
    skyvals[nvals] = sky[xi+1,yi]
    nvals += 1
    if yi ne 0 then begin
        skyvals[nvals] = sky[xi+1,yi-1]
        nvals+=1
    endif
    if yi ne maxY then begin
        skyvals[nvals] = sky[xi+1, yi+1]
        nvals += 1
    endif
endif
if yi ne 0 then begin
    skyvals[nvals] = sky[xi,yi-1]
    nvals += 1
endif
if yi ne maxY then begin
    skyvals[nvals] = sky[xi,yi+1]
    nvals += 1
endif
;print, nvals
;print, stddev(skyvals[0:nvals-1])
return, stddev(skyvals[0:nvals-1])

END


;NOTE: image here is made up of longs and has a bias of +1000, which
;needs to be subtracted to get the actual object counts
FUNCTION get_atlas_ivar, image, run, camcol, field, id, $
                         rerun=rerun, filter=filter, $
                         skyCount=skyCount, gain=gain,$
                         sky_err=sky_err, use_sky128=use_sky128

dimen=size(image, /dimensions )
ivar = dblarr(dimen[0], dimen[1])

fpObj = mrdfits( sdss_name('fpObjc', run, camcol, field, rerun=rerun), 1, rows=[id-1], /silent)
junk=mrdfits( sdss_name('fpObjc', run, camcol, field, rerun=rerun), 0, header, /silent)

;the dark varaince (read noise^2 + dark current is in the header
darkString = fxpar( header, 'DARK_VAR')
gainString = fxpar( header, 'GAIN' )

gain = -1.0D
dark_var = double((strsplit(darkString, /extract))[filter])
gain = double((strsplit(gainString, /extract))[filter])

;if something went wrong getting the gain:
if gain lt 1.0e-6 then begin
    print, 'could not find valid gain value, setting varaince to constant'
    ivar[where(image ne 1000)] = 1.0D
    return, ivar
endif

;the sky level associated iwth the object, this is the sky level
;calculated for that object, (has to do with the deblending properties
skyCount = fpObj.SKY[filter]
sky_err = fpObj.SKYERR[filter]

if (keyword_set(use_sky128)) then $
    sky_err = get_128sky_var(run, camcol, field, fpObj.ROWC[filter],$
                             fpObj.COLC[filter], $
                             rerun=rerun, filter=filter)

ivar = 1.0/(((image-1000)*1.0D + skyCount)/gain + dark_var+sky_err^2)
bad = where(ivar lt 1.0e-12, counts )
if counts ne 0 then ivar[where(ivar lt 1.0e-12)] = 0.0D
masked = where( image eq 1000 )
ivar[where(image eq 1000)] = 0.0D
return, ivar


END
