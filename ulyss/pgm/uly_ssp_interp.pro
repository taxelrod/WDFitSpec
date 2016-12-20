;+
; NAME:           ULY_SSP_INTERP
;
; PURPOSE:        Get interpolated spectrum for logT, Met and [Mg/Fe]
;
; USAGE:          model_new = uly_ssp_interp(model_grid, para)
;
; DESCRIPTION:    Interpolate spectrum from a grid of models at given
;                 log(age), [Fe/H] and [Mg/Fe]. Use ULY_2DERIV and
;                 INTERP_3D. Called, for example, from ULY_SSP_EXTR
;                 and ULY_FIT_LIN.
;
; ARGUMENTS:
;   model_grid:   Grid of model (structure)
;   [x, y, mgfe]: Log of age (Myr); metallicity (dex); [Mg/Fe]
;
; RETURN:         Interpolated spectrum (array)
;
; EXAMPLE:        Read a model grid and then extract a spectrum at
;                 given 10 Gyr and solar metallicity 
;                 model_grid=uly_ssp_read(uly_root+'/models/PHR_Elodie31.fits')
;                 model = uly_ssp_interp(model_grid, [alog(10000), 0.0])
;
; AUTHOR:         Philippe Prugniel
;-
; CATEGORY:       ULY_SSP
;-----------------------------------------------------------------------------

function uly_2deriv,x,y
;
; NAME:		ULY_2DERIV
;
; PURPOSE:      Compute the 2nd derivatives of Y(L,N) at points Xi (i=0,..N)
;
; USAGE:
;               y2 = uly_2deriv(x,y)
; DESCRIPTION:
;               Given arrays X of length N and Y of size LxN containing a
;               tabulated function, i.e.
;               Yi = f(Xi), with X1 < X2 < ... < Xn,
;               this routine returns and array Y2 of size LxN which contains
;               the second derivatives in the second coordinates of the
;               interpolating function at the tabulated points Xi.
;
;               Used by the functions ULY_SSP_INTERP and DERIVE_3D
; 
;
; ARGUMENTS:
;           x:  independent variable vector
;           y:  dependent variable vector
;
; RETURN:
;           y2: second derivatives at all x, of length n
;               the type of y2 is of the highest precision of x and y
;
; HISTORY:
;      D. Neill, October, 1991
;      http://www.astro.washington.edu/deutsch-bin/getpro/library43.html?SPLINT
;      adapted for ULySS,  Mina Koleva & Philippe Prugniel 20070308 
;-
;------------------------------------------------------------------------------

s = size(y)
n = s[2]

if s[2] ne n_elements(x) then message, 'X and Y should have the same length'

type = size(y, /TYPE) > size(x, /TYPE)
y2 = make_array(s[1],n, /NOZERO, TYPE=type)
u = make_array(s[1],n, /NOZERO, TYPE=type)
;
; The lower boundary condition is set to be "natural"
;
y2[*,0] = 0.
u[*,0] = 0.

;
; This is the decomposition loop of the tridiagonal algorithm.  Y2 and
; U are used for temporary storage of the decomposed factors.
; We make the choice of using loops along the interpolating region
; rather than using TRISOL or SPL_INIT, because this allows to use 
; array operations along the other axis. This is the best solution because
; this other axis is in principle much longer (i.e. it is the number
; of wavelength points and is >10 times larger than the length of the
; interpolated axis).
;
tsig = double((x - shift(x, -1))) / (shift(x, +1) - shift(x, -1))

status = check_math(MASK=32)  ; this function may generate underflows

for i=1,n-2 do begin
    p = tsig[i] * y2[*,i-1] + 2.
    y2[*,i] = (tsig[i]-1.) / p

    u[*,i]=( 6. * ((y[*,i+1]-y[*,i]) / (x[i+1]-x[i]) - (y[*,i]-y[*,i-1]) $
                   / (x[i]-x[i-1])) / (x[i+1]-x[i-1]) - tsig[i]*u[*,i-1]) / p
endfor

; The upper boundary condition is set either to be "natural"
;
y2[*,n-1] = 0.
;
; This is the backsubstitution loop of the tridiagonal algorithm
;
for k=n-2,0,-1 do begin
    y2[*,k] = y2[*,k] * y2[*,k+1] + u[*,k]
endfor

if status eq 0 then status = check_math(MASK=32)  ; dismiss underflows

return, y2
end


;-----------------------------------------------------------------------------
function interp_3d, x, y, TMPARR=tmparr, D2X=d2x,$
                    XGRID=xgrid, YGRID=ygrid
compile_opt idl2, hidden

; NAME:         INTERP_3D (internal routine)
;
; PURPOSE:      Interpolate in a data cube
;
; USAGE:        z = interp_3d( x, y, TMPARR=tmparr, D2X=d2x,$
;                               XGRID=xgrid, YGRID=ygrid)
; DESCRIPTION:  Interpolate in a data cube (tmparr) at (x,y) coordinates
;
; ARGUMENTS:
;      x, y:    The demanded coordinates at which the cube should be
;               intepolated
;
; KEYWORDS:
;   TMPARR :    model cube 
;   D2X    :    array of X 2nd derivatives (input)
;   XGRID  :    array of X for the cube (input)
;   YGRID  :    array of Y for the cube (input)
;
; HISTORY:       <unknown> author

if not finite(x) then message, 'Input (X) is not finite'

s = size(tmparr)

if(x gt max(xgrid)) then x = max(xgrid)
if(y gt max(ygrid)) then y = max(ygrid)
if(x lt min(xgrid)) then x = min(xgrid)
if(y lt min(ygrid)) then y = min(ygrid)

if(s[0] eq 1) then return, tmparr

;    t=systime(1)

if((s[0] eq 3) and keyword_set(d2x) and (s[3] ge 4)) then begin

; determine in what region is the point[x,y] 
    xval = interpol(findgen(s[2]), xgrid, x) ;
    yval = interpol(findgen(s[3]), ygrid, y)
    xmin = floor(xval) 
    xmax = xmin + 1
    if(xmax ge s[2]) then xmax=xmin
    ymin = floor(yval)
    ymax = ymin + 1
    if(ymax ge s[3]) then ymax=ymin
    ys1 = ymin - 1              ;limit for the spline
    if(ys1 lt 0) then ys1 = ys1 + 1
    ys2 = ys1 + 3    ;we need 4 points of metallicities to make spline
    if(ys2 ge s[3]) then begin
        ys2 = s[3] - 1
        ys1 = ys2 - 3
    endif
    
    intvec = dblarr(s[1], /NOZERO) ;init. interpolated vector

;interpolate along x direction for the 4 y nodes that we need for the
;spline in y

    h = (xgrid[xmax]-xgrid[xmin])
    ax = (h eq 0d) ? 0 : (xgrid[xmax]-x) / h
;    ax = (xgrid[xmax]-x) / (xgrid[xmax]-xgrid[xmin])
    bx = 1.0d - ax
    cx = (ax^3-ax) * ((xgrid[xmax]-xgrid[xmin])^2) / 6d
    dx = (bx^3-bx) * ((xgrid[xmax]-xgrid[xmin])^2) / 6d
    mm = reform(ax*tmparr[*,xmin,ys1:ys1+3]+bx*tmparr[*,xmax,ys1:ys1+3]+$
                cx*d2x[*,xmin,ys1:ys1+3]+dx*d2x[*,xmax,ys1:ys1+3], s[1], 4)
    d2y = uly_2deriv(ygrid[ys1:ys1+3],mm)

    h = (ygrid[ymax]-ygrid[ymin])
    ay = (h eq 0) ? 0 : (ygrid[ymax]-y)/h
;    ay = (ygrid[ymax]-y)/h
    by = 1.0d - ay
    cy = ay^3 - ay
    dy = by^3 - by
    intvec[*] = ay*mm[*,ymin-ys1]+by*mm[*,ymax-ys1]+$
          (cy*d2y[*,ymin-ys1]+dy*d2y[*,ymax-ys1])*(1d/6d)*(h^2)
endif else begin
    if s[0] ne 3 then message, 'The model grid must be 3D'
    if s[3] le 4 then message, 'The grid must have > 3 points for spline interpolation'
endelse

; print, 'TIME INTERP 3D',systime(1)-t

return, intvec
end

; ----------------------------------------------------------------------------
function uly_ssp_interp, model_grid, para
compile_opt idl2, hidden; para[0] : log age
; para[1] : metallicity
; para[2] : Mg/Fe

; Diagnostic required extrapolations
m0 = min(*model_grid.o_age, MAX=m1)
if para[0] lt m0 or para[0] gt m1 then message, /INFO, 'Requested age is out of range ... extrapolation'
m0 = min(*model_grid.o_metal, MAX=m1)
if para[1] lt m0 or para[1] gt m1 then message, /INFO, 'Requested [Fe/H]='+$
  strtrim(para[1],2)+' is out of the valid range ('+strtrim(m0,2)+$
  ','+strtrim(m1,2)+')  ... extrapolation'

if size(*model_grid.data, /N_DIM) gt 3 then begin  ; grid with Mg/Fe resol
    if size(*model_grid.data, /N_DIM) ne 4 $
      or (size(*model_grid.data, /DIM))[3] ne 2 then $
      message, 'Support on grid with 2 [Mg/Fe] levels)'
    inpmgfe = *model_grid.o_mgfe
    const = (exp(inpmgfe[1])-exp(para[2])) / (exp(inpmgfe[1])-exp(inpmgfe[0]))
    return, const * $
      interp_3d(para[0], para[1], TMPARR=(*model_grid.data)[*,*,*,0], $
                D2X=(*model_grid.d2t)[*,*,*,0],$
                XGRID=*model_grid.o_age, YGRID=*model_grid.o_metal) + $
      (1. - const) * $
      interp_3d(para[0], para[1], TMPARR=(*model_grid.data)[*,*,*,1], $
                D2X=(*model_grid.d2t)[*,*,*,1],$
                XGRID=*model_grid.o_age, YGRID=*model_grid.o_metal)
    
endif else $
  return, interp_3d(para[0], para[1], TMPARR=*model_grid.data, $
                    D2X=*model_grid.d2t,$
                    XGRID=*model_grid.o_age, YGRID=*model_grid.o_metal)
end

; ---- end -------------------------------------------------------------------
