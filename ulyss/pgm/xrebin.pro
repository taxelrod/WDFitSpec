;+
; NAME:              XREBIN
;
; PURPOSE:           Rebin a vector
;
; USAGE:
;         yout = xrebin(xin, yin, xout [,/SPLINE] [,/LSQUADRATIC] 
;                       [,/QUADRATIC], [/SPLINF], [/CUBIC], [/NAN]
; ARGUMENTS:
;         XIN         (N+1) points of the borders of each input bin
;         YIN         values (N) of input bins (1 or 2D)
;         XOUT        points of the borders of output bins
; KEYWORDS:
;         SPLINE      Spline interpolation (default: linear)
;         LSQUADRATIC Least square quadratic interpolation
;         QUADRATIC   Quadratic interpolation
;         SPLINF      Other version of spline interpolation
;                     3-4 times faster than SPLINE, the results are
;                     slightly different, but seemingly better.
;         CUBIC       Use the function INTERPOLATE with CUBIC=-05
;                     Cubic convolution is an interpolation method that  
;                     closely approximates the theoretically optimum sinc  
;                     interpolation function using cubic polynomials.
;         NAN         To treat the missing values in YIN as
;                     zeros when interpolating it. Use with caution. 
; OUTPUT
;         YOUT        output 1 or 2D array, rebined at the borders
;                     requested by XOUT. 
;                     YOUT is in double precision (it shall be converted
;                     in float if appropriate.
;
; DESCRIPTION:   
;    This function rebins a vector or a n-dimensional array into
;    different bins along the first axis.
;
;    The rebining is made as integrating the signal, interpolating and
;    differentiating the result. The integrated signal is preserved.
;
;    Different interpolating algorithms are proposed, and only one of
;    the keywords SPLINE LSQUADRATIC QUADRATIC SPLINF CUBIC must be given.
;    If not of them is specified, the routine INTERPOL will make a linear
;    interpolation. This latter routin also handle the SPLINE LSQUADRATIC
;    QUADRATIC cases. The difference between SPLINF and SPLINE is
;    that the spline representation is computed on the whole vector
;    instead of on 4-points blocks (this is faster and often our preferred
;    option). CUBIC is using the function INTERPOLATE with CUBIC=-0.5 (sinc).
;
;    No check is performed to prevent extrapolation (beware to give proper input).
;
;    The computations are systematically made in double precision. The
;    reason is that the algorithm results in loss of relative precision equal
;    to the number of pixels. In float, this is often a unacceptable.
;
; AUTHOR: Mina Koleva (2007/03/01)
;-
; CATEGORY:    ULY_UTIL
;------------------------------------------------------------------------------
function xrebin, xin, yin, xout, SPLINE=spline, LSQUADRATIC=lsquadratic, QUADRATIC=quadratic, SPLINF=splinf, CUBIC=cubic, NAN=nan

compile_opt idl2
on_error, 2

s = size(yin)
if s[0] lt 1  then message,'Yin should be a vector or array'

if s[1]+1 ne (size(xin))[1] then message,'Xin must have one element more than the 1st dim of Yin '

if s[0] eq 1 then begin    ; case of a 1D array

    integr = [0, total(double(yin), 1, /CUMULATIVE, NaN=NaN)]  ; integrate yin
;   Note that we force the computation in double precision
;   The integration-differentiation process results in a loss of precision
;   which is often critical in float (and acceptable in double)
;   The loss of precision is equal to the number of pixels, so for 1e6
;   pixels, if the computation were made in float, the resulting precision
;   would not be better than a few 10%...

    ;interpolate yout for xout
    case (1) of
        keyword_set(splinf): begin
;           if spl_init generates an underflow error we clear it
            status = check_math(MASK=32)
            y2 = spl_init(xin, integr)
            integr_interp = spl_interp(xin, integr, y2, xout)
;           If any of the 4 input vectors is double precision, the
;           result is also in double precision
            if status eq 0 then status = check_math(MASK=32)
        end

        keyword_set(cubic): begin
            iout = interpol(lindgen(n_elements(xin)), xin, xout)
            integr_interp = interpolate(integr, iout, /GRID, CUBIC=-0.5)
        end

        else : $
          integr_interp = $
          interpol(integr, xin, xout, $
                   SPLINE=spline, LSQUADRATIC=lsquadratic, QUADRATIC=quadratic)
        
    endcase
    ;differentiate
    return, (shift(integr_interp,-1)-integr_interp)[0:n_elements(integr_interp)-2]
    
endif else begin   ; case of a nD array
    dim = size(yin, /DIM)
    integr = [dblarr([1,dim[1:*]]), total(double(yin), 1, /CUMULATIVE, NAN=NaN)]  ; integrate yin
    integr = reform(integr, dim[0]+1, n_elements(integr)/(dim[0]+1), /OVER)
    sr = size(integr)

    case (1) of
        keyword_set(splinf): begin
            integr_interp = make_array(n_elements(xout), sr[2], /DOUB, /NOZERO)
            status = check_math(MASK=32)
            for i=0, sr[2]-1 do begin
                y2 = spl_init(xin, integr[*,i])
                integr_interp[*,i] = spl_interp(xin, integr[*,i], y2, xout)
            endfor
            if status eq 0 then status = check_math(MASK=32)
        end

        keyword_set(cubic): begin
            iout = interpol(lindgen(n_elements(xin)), xin, xout)
            yout = lindgen(sr[2])
            integr_interp = interpolate(integr, iout, yout, /GR, CUB=-0.5)
        end

        else: begin
; note: we can certainly improve this significantly if we avoid the
;  VALUE_LOCATE in interpol
            integr_interp = make_array(n_elements(xout), sr[2], /DOUB, /NOZERO)
            for i=0, sr[2]-1 do $
              integr_interp[*,i] = $
              interpol(integr[*,i], xin, xout, $
                       SPLINE=spline, LSQUADRATIC=lsquadratic, QUAD=quadratic)
        end
    endcase
    ;differentiate
    s1=(size(integr_interp))[1]
    return, reform((shift(integr_interp,-1)-integr_interp)[0:s1-2, *], $
                   [s1-1, dim[1:*]])
endelse

; note: We can also use the function INTERPOLATE (with /GRID,CUB=-0.5)
;  We would not need the for loop (in the 2D case)
;  

end

;------- TESTING XREBIN----------
pro test_xrebin

;test 1
print, 'TEST 1...',FORMAT="(A,$)"
xin = [1,2,3,4] 
yin = [1,2,1]
xout= [1,2,3,4]
yout = xrebin(xin,yin,xout)
a = where(yout - [1,2,1],n)
if n eq 0 then $
print, 'OK' else $
print, 'FAILED'

;test 2
print, 'TEST 2...',FORMAT="(A,$)"
xin = [1,2,3,4] 
yin = [1.,2.,1.]
xout = [1,1.5,2,2.5,3,3.5,4]
yout = xrebin(xin,yin,xout)
a = where(yout - [.5,.5,1,1,.5,.5],n)
if n eq 0 then $
print, 'OK' else $
print, 'FAILED'

;test 3
print, 'TEST 3...',FORMAT="(A,$)"
xin = [1,1.5,2,2.5,3,3.5,4]
yin = [0.5,0.5,1.,1.,.5,.5] 
xout = [1,2,3,4]
yout = xrebin(xin,yin,xout)
a = where(yout - [1,2,1],n)
if n eq 0 then $
print, 'OK' else $
print, 'FAILED'


end

;--- end ---

