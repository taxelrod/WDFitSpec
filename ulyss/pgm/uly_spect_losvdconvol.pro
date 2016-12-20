;+
; NAME:      
;                 ULY_SPECT_LOSVDCONVOL
;
; PURPOSE:   
;                 Convolve a spectrum with a losvd
;
; USAGE:
;                 spectr = uly_spect_losvdconvol(SignalIn, cz, sigma, h3, h4, /OVERWRITE)
;
; ARGUMENTS:
;   SignalIn:     Input data structure (spectrum)
;
;   cz:           Input systemic cz in km/s
;
;   sigma:        Input velocity dispersion in km/s
;
;   h3:           Input 3rd Hermite coefficient
;
;   h4:           Input 4rd Hermite coefficient
;
; KEYWORDS:
;   OVERWRITE     Set this keyword if input and output spect are the
;                 same variable as to prevent a memory leak
;
; OUTPUT:
;                 Spectrum structure containing the convolved data
;
; DESCRIPTION:
;     The input data must be a 'spectrum structure', and the data
;     can be sampled linearly or in log (see ULY_SPECT_ALLOC).
;     If the sampling is not logarithmic, the routine first rebin
;     the input before convolving, and the output spectrum is log sampled.
;
;     The error spectrum and pixel mask, if present in (SignalIn) are
;     also transformed in accordance with the signal. 
;
;     If the data array contained in (SignalIn) has two dimensions,
;     each line is convolved similarly. 
;
; HISTORY:        Ph. Prugniel (2008/01/31) created
;-
; CATEGORY:       ULY
;------------------------------------------------------------------------------
; Hermite Polynoms   (from Cappellari)
; see Marel and Franx, 1993, ApJ, 407
; (Annexe A)
;------------------------------------------------------------------------------
function Hermite3, x
h = (2*x^3 - 3*x)/sqrt(3.)
return, h
end

function Hermite4, x
h = (4*x^4 - 12*x^2 + 3)/sqrt(24.)
return, h
end

;------------------------------------------------------------------------
; From Cappellari, bug fix Prugniel (underflow)
; Compute the LOSVD profile which have to be convolved with the spectrum
; see Cappellari and Emsellem, 2004, PASP, 116
;
; ARGUMENTS :
;   cz:           Systemic cz in km/s
;
;   sigma:        Dispersion velocity in km/s. Defaulted to 0.
;
;   h3, h4:       Gauss-Hermite coefficients. Defaulted to 0.
;
;   resol:        Size (in km/s) of each pixel of the output array
;
; OUTPUTS :
;   slitout:      Output array
;
pro uly_slit, cz, sigma, h3, h4, velscale, slitout

c = 299792.458d
logl_shift = alog(1d + cz/c) / velscale * c     ; shift in pixels
logl_sigma = alog(1d + sigma/c) / velscale * c  ; sigma in pixels

N = ceil((abs(logl_shift) + 5d*logl_sigma))
x = N - dindgen(2*n+1)
y = (x-logl_shift)/logl_sigma
large = where(abs(y) gt 30, c)
if c gt 0 then y[large] = 30d  ; prevent underflow error
slitout = EXP(-y*y/2d) / logl_sigma / sqrt(2d*!pi) 
slitout = slitout*( 1d + h3*Hermite3(y) + h4*Hermite4(y))
slitout /= total(slitout)

end

function uly_spect_losvdconvol, SignalIn, cz, sigma, h3, h4, OVERWRITE=overwrite

compile_opt idl2
on_error, 2

; first rebin in log if necessary.
input_sampling = SignalIn.sampling
if input_sampling ne 1 then SignalTmp = uly_spect_logrebin(SignalIn) $
else SignalTmp = SignalIn

if n_elements(cz) eq 0 then $
   message, 'ULY_SPECT_LOSVDCONVOL: Velocity (cz) must be given'
if n_elements(sigma) eq 0 then sigma = 0
if n_elements(h3) eq 0 then h3 = 0
if n_elements(h4) eq 0 then h4 = 0

s = size(*SignalTmp.data)
if s[0] ne 1 then $
  message, 'This function process only 1D arrays (yet)'


c = 299792.458d                 ; Speed of light in km/s

velscale = SignalTmp.step * c

if keyword_set(overwrite) then SignalOut=SignalTmp else $
   SignalOut = uly_spect_alloc(SPECTRUM=SignalTmp)

uly_slit, cz, sigma, h3, h4, velscale, losvd

status = check_math(MASK=32, /NOCLEAR)
*SignalOut.data = convol((*SignalTmp.data),losvd,/EDGE_TRUNCATE)
if status eq 0 then status = check_math(MASK=32)

if n_elements(*SignalTmp.err) ne 0 then begin
    if n_elements(*SignalTmp.err) ne n_elements(*SignalTmp.data) then $
      message, 'spectrum structure unconsistent (err)'
    status = check_math(MASK=32, /NOCLEAR)
    err = convol((*SignalTmp.err)^2, losvd, /EDGE_TRUNCATE)
    if status eq 0 then status = check_math(MASK=32)
    ; the smoothing reduces the variance by the factor 1/total(losvd^2)
    *SignalOut.err = sqrt(err * total(losvd^2))
endif

if n_elements(*SignalTmp.goodpix) ne 0 then begin
; need to (1) make a spectrum for the pixel list, (2) convolve it (3)
; remake a list
    maskI = replicate(0E, n_elements(*SignalTmp.data))
    maskI[*SignalTmp.goodpix] = 1
    maskI = reform(maskI, size(*SignalIn.data, /DIM), /OVER)

;   If the convolution generates an underflow, we clear the signal
    status = check_math(MASK=32, /NOCLEAR)
    maskO = convol(maskI, losvd, /EDGE_TRUNCATE)
    if status eq 0 then status = check_math(MASK=32)
    *SignalOut.goodpix = where(maskO gt 0.5, cnt)
    if cnt eq 0 then message, $
      /INFO, 'No good pixels were left after the convolution '+ $
      '(it will probably be interpreted as "all pixels are good")'
endif

; The degree-of-freedom factor has to be changed ... the number of
; independent' points decreases.
dof_factor = 1./max(losvd)
; We may think of improving the computation of dof_factor...
; dof_factor *= 1.5 ; dof_factor is too small for large sigmas with this method 

SignalOut.dof_factor = SignalTmp.dof_factor * dof_factor

sxaddpar, *SignalOut.hdr, 'HISTORY', 'uly_spect_losvdconvol'

input_sampling = SignalIn.sampling
if input_sampling ne 1 then SignalTmp = uly_spect_logrebin(SignalIn) $
else SignalTmp = SignalIn

if input_sampling ne 1 then begin
   if  keyword_set(overwrite) then uly_spect_free, SignalIn $
   else uly_spect_free, SignalTmp
endif

return, SignalOut
end
