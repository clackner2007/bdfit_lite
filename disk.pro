;-----------------------------
;  Claire Lackner
;  Jan 17 2010
;
;  forms an image of a disk galaxy with a given scale length, scale
;  height, inclination, intensity, and integration cutoff radius
;
;  the images are also recentered and rotated, all of these parameters
;  are stored in the 1x8 array params
;
;  a_norm      normalization of surface brightness profile ==
;               4*pi*r0^2*H0*rho0 / (pi*R0^2) = 4*rho0*H0
;  R0          disk scale length
;  H0          disk scale height
;  Inc         Inclination of disk
;  cutRad      radius at which integration stops, z direction
;  x_0         central position in x
;  y_0         central position in y
;  phi         rotation angle of major axis
;
;
;-----------------------------


;used in the the density functions, rho and rhoarray
;returns the inclined y-value
FUNCTION yinc, y, z, i
yinc = y*cos(i) + z*sin(i)
return, yinc
END

;returns inclided z-value
FUNCTION zinc, y, z, i
zinc = -1.*y*sin(i) + z*cos(i)
return, zinc
END

;-----------------------------
;-----------------------------
;functions returns the 3D density that needs to be integrated over 
FUNCTION rho, z, _EXTRA=_EXTRA

rho0 = _EXTRA.norm
R0 = _EXTRA.r0
H0 = _EXTRA.h0
i = _EXTRA.inc
x = _EXTRA.x
y = _EXTRA.y

y = dblarr( n_elements(z) ) + _EXTRA.y
x = dblarr( n_elements(z) ) + _EXTRA.x


yy=yinc(y,z,i)

rho = rho0*exp(-(x*x+yy*yy)^(0.50)/R0 - abs(zinc(y,z,i))/H0)

return, rho

END

;functions returns the 3D density that needs to be integrated over
;assumes that x and y are arrays, and that qsimparray is called
FUNCTION rhoarray, z, _EXTRA=_EXTRA

rho0 = _EXTRA.norm
R0 = _EXTRA.r0
H0 = _EXTRA.h0
i = _EXTRA.inc
num = n_elements(z)
;put x and y into a set of 3-d matrices, one array for each entry in z
;uses cmreplicate in ~/pro, generalized replicate
x = cmreplicate( _EXTRA.x, [num] )
y = cmreplicate( _EXTRA.y, [num] )
nx = (size(x, /dimensions))[0]
ny = (size(x, /dimensions))[1]

if (size(z))[0] eq 0 then zmat = cmreplicate( [z], [nx, ny] ) else $
  zmat = cmreplicate( z, [nx,ny] )
zmat = transpose(zmat,[1,2,0])
yy=yinc(y,zmat,i)

rho = rho0*exp(-(x*x+yy*yy)^(0.50)/R0 - abs(zinc(y,zmat,i))/H0)

return, rho

END

;-----------------------------
;-----------------------------

;returns surface brightness image for a double exponential disk galaxy
FUNCTION disk, x, y, params

x0 = params[5]
y0 = params[6]
phi = params[7]
inc = params[3]

;translate + rotate
xtrans = (x - x0) * cos(phi) + (y - y0) * sin(phi)
ytrans = (y - y0) * cos(phi) - (x - x0) * sin(phi)


;do each integral one by one, for now, it turns out to be faster, the
;amount of computation saved by not doing things past the tolerance is
;bigger than the amount killing one for loop buys you
  image = dblarr( (size(x,/dimensions))[0], (size(x,/dimensions))[1] )

  for i=0L, n_elements(xtrans)-1 do begin
      extras = {norm:params[0]/(4.0*params[2]), r0:params[1], $
                h0:params[2], inc:inc, x:xtrans[i], y:ytrans[i]}

      ;guess of integral = 2*Rcut*rho(z=0)
      guess=2.0*params[4]*rho(0.0, _EXTRA=extras )
    
      qsimp, 'rho', -1.0*params[4], 1.0*params[4], guess, $
        _EXTRA=extras

      image[i] = guess
    
  endfor
;does everything in arrays
;  extras = {norm:params[0]/(4.0*params[2]), r0:params[1], $
;            h0:params[2], inc:inc, x:xtrans, y:ytrans}
 
;  image = 2.0*params[4]*rhoarray(0.0, _EXTRA=extras )

;  qsimparray, 'rhoarray', -1.0*params[4], 1.0*params[4], image, _EXTRA=extras


return, image

END

