;+
; NAME:      
;                 ULY_SPECT_LOGREBIN
;
; PURPOSE:  
;                 Rebin logarithmically a spectrum 
;
; USAGE:
;                 slog = uly_spect_logrebin(<spectrum> [,<vsc>])
;
; ARGUMENTS:
;   <spectrum>:   Input, data structure (spectrum).
;
;   <vsc>:        I/O, size of each pixel in km/s for the output data.
;
; KEYWORDS:
;   WAVERANGE     Input scalar or 2-elements vector. Wavelength range.
;                 If a single value s passed, it is interpreted as the
;                 start of the wavelength range.
;
;   FLUX_CONSERVE Input boolean.
;                 Set this keyword to skip the intensity correction
;                 for the change of the local pixel size.
;                 If flux is conserved, intensity is not; by default
;                 intensity (i.e. flux density) is conserved.
;
;   LINEAR        Input boolean
;                 Set this keyword to use a linear interpolation of the
;                 cumul function. By default a spline interpolation is made.
;                 A linear interpolation implies a convolution by a one-pixel
;                 width boxcar function, while a spline interpolation
;                 performs an implicite deconvolution of the pixels signature.
;                 The difference may be sentitive in case of rebinning to
;                 similar or finer pixels.
;
;   EXACT         Input boolean.
;                 Specify that the output must match exactly the WCS
;                 specified by <waverange> and <vsc>. Otherwise approximations
;                 can be made to avoid shifts by small fractions of pixels.
;
;   OVERWRITE     Input boolean.
;                 Set this keyword if input and output spect are the
;                 same variable, to save memory and prevent a memory leak.
; OUTPUT:
;                 SLOG, spectrum structure containing the log rebinned data.
;
; DESCRIPTION:
;     This routine performs rigorously a rebinning into logarithmically
;     spaced wavelength points. The algorithm consists in:
;     (i) integrate the spectrum, (ii) resample (i. e. interpolate)
;     this integral and (iii) differentiate the resampled integral. This
;     naturally conserves the flux. The resampling is done with a spline
;     or linear interpolation (see XREBIN).
;
;     As it conserves the flux, the intensity of each pixel change
;     because the size of a pixel change. As well, the shape of the
;     intensity distribution changes, for example when converting a 
;     linearly sampled spectrum into a log-sampled one. To alleviate
;     this aspect, the function corrects by default for the change of  
;     the pixel size. in is way the unit of the Intensity-axis are preserved.
;     To skip this correction, use the keyword FLUX_CONSERVE.
;
;     The input data must be a 'spectrum structure' (see ULY_SPECT_ALLOC), 
;     and the data can be sampled linearly, in log or at explicit wavelength.
;     
;     By default the routine avoids a resampling if the input spectrum
;     has the required step (<vsc>), and in that case simply extract
;     the spectral range specified with <waverange>. Therefore, in this
;     case the spectrum may be shifted by a fraction of pixel with respect
;     to the requirement. To force a resampling to the exact grid,
;     specify the keyword /EXACT. (note that in the case the constraint on the
;     wavelength of the last pixel cannot be exact).
;
;     <vsc> is the velocity scale, ie. size in km/s of each pixel in the
;     log-rebinned spectrum.
;     If <vsc> is initialized in input, this value will be used to rebin
;     the input spectrum, otherwise the velocity scale is computed
;     in order to keep the same number of pixels in output as in input.
;     The computed value is returned in <vsc>.
;
;     WAVERANGE. If a scalar or 1-element array is passed, it is the
;     wavelength of the center of the output spectrum's 1st pixel in Angstrom. 
;     If a 2-elements array array is passed, the wavelength range is [start, end]
;     If the specified start is greater than the last wavelength of the 
;     input spectrum, the routine fails with an error. 
;     If the specified start is smaller than the begining of the spectrum,
;     the effective start is shifted by a round number of pixels.
;     If <waverange> is not passed, the spectrum is rebinned with a WCS
;     so that the beginning edge of the first pixel coincide with each
;     other both in input and output spectra.
;
;     The error spectrum and pixel mask, if present in (spectrum) are
;     also transformed in accordance with the signal. 
;
;     If the data and error arrays contained in (spectrum) have two dimensions,
;     all the lines are rebinned similarly. 
;
; SEE ALSO:       ULY_SPECT_LINREBIN, ULY_SPECT_READ
;
; REQUIRED ROUTINES:  
;                 XREBIN.PRO
;
; HISTORY:        Mina Koleva (2007/03/01)
;     
;-
; CATEGORY:       ULY
;------------------------------------------------------------------------------
function uly_spect_logrebin, SignalIn, vsc, WAVERANGE=waverange, $
                             FLUX_CONSERVE=flux_conserve,        $
                             LINEAR=linear, EXACT=exact, OVERWRITE=overwrite

compile_opt idl2
on_error, 2

if size(SignalIn, /TYPE) ne 8 then message, 'Input must be a spectrum structure'
if n_elements(*SignalIn.data) eq 0 then message, 'Input spectrum undefined'

npix = (size(*SignalIn.data,/DIM))

c = 299792.458d                 ; Speed of light in km/s

if keyword_set(linear) then splinf = 0 else splinf = 1

if SignalIn.sampling eq 1 then begin
    cpout = 0
;   Note that comparison between SignalIn.step and vsc/c tolerates
;   rounding errors (0.001 px over the whole range)
    if (n_elements(vsc) eq 1) then begin
        if (abs(1d - double(vsc)/c/SignalIn.step)*npix[0]) lt 0.001 then cpout = 1
    endif else if n_elements(vsc) eq 0 then cpout = 1
    if keyword_set(exact) and cpout eq 1 and n_elements(waverange) gt 0 then begin
        nshift = (SignalIn.start-alog(waverange[0]))/SignalIn.step mod 1d
        if abs(nshift) gt 0.001 then cpout = 0

        if cpout eq 1 and n_elements(waverange) eq 2 then begin
           nshift = (SignalIn.start + (npix[0]-1)*SignalIn.step $
                    - alog(waverange[0]))/SignalIn.step mod 1d
           if abs(nshift) gt 0.001 then cpout = 0
        endif
    endif
    if cpout eq 1 then begin    ; just have to copy input to output...
        vsc = c * SignalIn.step ; for when vsc arg is not set
        return, uly_spect_extract(SignalIn, WAVERANGE=waverange, $
                                  OVERWRITE=overwrite)
    endif
endif

;compute the borders of input bins
if SignalIn.sampling eq 0 then begin
    indxarr = lindgen(npix[0])
    borders = SignalIn.start + [-0.5d,$
                                (indxarr+0.5d)] * SignalIn.step
    bordersLog = alog(borders)

endif else if SignalIn.sampling eq 1 then begin
    indxarr = lindgen(npix[0])
    bordersLog = SignalIn.start + [-0.5d,$
                                   (indxarr+0.5d)] * SignalIn.step
    borders = exp(bordersLog)
endif else if SignalIn.sampling eq 2 then begin
    borders = [(*SignalIn.wavelen + shift(*SignalIn.wavelen, 1)) / 2d, 0d]
    borders[0] = 2.*(*SignalIn.wavelen)[0] - borders[1]
    borders[npix[0]] = 2.*(*SignalIn.wavelen)[npix[0]-1] - borders[npix[0]-1]
    bordersLog = alog(borders)
endif else $
  message, 'Sampling mode of input is invalid:'+SignalIn.sampling

; Determine the velocity scale, vsc, preserving the total number of pixels
if n_elements(vsc) eq 0 then begin
    if SignalIn.sampling eq 1 then vsc = SignalIn.step * c $
    else if SignalIn.sampling eq 2 then begin
;       We assume that the wavelength are ordered, and we do not check it
        if n_elements(waverange) gt 0 then begin
            nw = value_locate(*SignalIn.wavelen, waverange)
            if nw[0] eq -1 then nw[0] = 0
            if n_elements(waverange) eq 1 then nw = [nw[0], npix[0]]
            if nw[1] eq -1 then nw[1] = 0 ; this is bad
        endif else nw = [0d, npix[0]-1]
        wr = (*SignalIn.wavelen)[nw]
        vsc =  alog(wr[1]/wr[0]) / npix[0] * c 
    endif else begin
        if n_elements(waverange) gt 0 then begin
            nw = (waverange - SignalIn.start)/SignalIn.step -0.5d
            if n_elements(waverange) gt 1 then $
              nw[1] = min([nw[1]+1, npix[0]-0.5]) $
            else nw = [nw[0], npix[0]-0.5]
        endif else nw = [-0.5d, npix[0]-0.5]

        wr = alog(SignalIn.start + nw * SignalIn.step)
        vsc =  (wr[1]-wr[0]) / (nw[1] - nw[0]) * c 
    endelse
endif

logScale = vsc/c   ; Scaling in the output spectrum

; Determine the start wavelength, logStart, and number of pix, nout, in output
logRange = SignalIn.start + [-0.5d,(npix[0]-0.5)] * SignalIn.step
if SignalIn.sampling eq 0 then logRange = alog(logRange) else $
if SignalIn.sampling eq 2 then logRange = alog([(*SignalIn.wavelen)[0],(*SignalIn.wavelen)[npix-1]])
logRange += [0.5d,-0.5d]*logScale
if n_elements(waverange) gt 0 then begin
;  the 1D-7 is a tolerance for rounding errors
   nshift = ceil(max([0d, (logRange[0] - alog(waverange[0])) / logScale - 1d-7]))
   logRange[0] = alog(waverange[0]) + logScale * nshift  ; center of 1st out pix

   if n_elements(waverange) eq 2 then $
      logRange[1] = min([alog(waverange[1]), logRange[1]]) 

   if logRange[1] lt logRange[0] then begin
       message, /INFO, 'waverange is not valid: '+ $
         strjoin(strtrim(waverange,2),',')
       return, 0
   endif
endif 
nout = round((logRange[1]-logRange[0]) / logScale + 1)
logStart = logRange[0]

; determine the DOF factor
; Note that this is an approximation. Even if the number of pixels is conserved
; the rebinning introduces a correlation. We shall improve this...
if SignalIn.sampling lt 2 then begin
   if SignalIn.sampling eq 0 then logRange = exp(logRange)
   nin = round((logRange[1]-logRange[0]) / SignalIn.step + 1)
endif else nin = n_elements(*SignalIn.wavelen)
dof_factor = nout/double(nin)

; we may think of a correction of dof_factor in case of splinf (obvious reason) 
; (the following one is OK only for large factors, it is too strong for small factors:)
;if splinf eq 1 and dof_factor gt 1 then dof_factor *= 0.86

if logStart-logScale/2d gt bordersLog[npix[0]] then begin
    message, 'start value is not in the valid range', /INFO
    return, 0
endif

; define the new borders of the bins
NewBorders = exp(logStart + (dindgen(nout+1)-0.5d) * logScale)

dim = size(*SignalIn.data, /DIM)
n_data = n_elements(*SignalIn.data)
n_err = n_elements(*SignalIn.err)
n_dim = size(*SignalIn.data, /N_DIM)

statmath = check_math(MASK=32, /NOCLEAR)  ; do we have an underflow error on the stack?

; Determine the conversion factor/vector accounting for pixel size change
case SignalIn.sampling of
    0: begin
        flat = exp(logStart + dindgen(nout)*logScale) * logScale/SignalIn.step
        if n_dim gt 1 then flat = rebin(flat, [nout, dim[1:*]])
    end
    1: begin        
        flat = logScale/SignalIn.step
    end
    else: begin
       flat = fltarr(dim[0], /NOZERO) ; faster to create an array and
       replicate_inplace, flat, 1B    ; populate with 1s 
       flat = xrebin(borders, flat, NewBorders, SPLINF=splinf) 
       if n_dim gt 1 then flat = rebin(flat, [nout, dim[1:*]])
    end
endcase

if keyword_set(overwrite) then SignalOut = SignalIn $
else SignalOut = uly_spect_alloc()

SignalIn_type = size(*SignalIn.data, /TYPE)

if keyword_set(flux_conserve) then $
  *SignalOut.data = xrebin(borders,*SignalIn.data,NewBorders,SPLINF=splinf) $
else $
  *SignalOut.data = xrebin(borders,*SignalIn.data,NewBorders,SPLINF=splinf) / flat

if SignalIn_type lt 5 then *SignalOut.data = float(temporary(*SignalOut.data))

if n_err eq n_data then begin
; if err is very small, then the ^2 is making it close to 0 and log_rebinning
; may produce negative values. Even small values can have dramatic 
; consequences, as they give an exagerated weight to the affected pixels.
; To prevent subsequent errors we set these small values to
; the min error in the input spectrum.
; It may possibly be better to simply use a linear interpolation
; that would not induce such artifacts...
    err = xrebin(borders,*SignalIn.err^2,NewBorders,SPLINF=splinf) 
    if n_elements(*SignalIn.goodpix) ne 0 then $
       minerr = min((*SignalIn.err)[*SignalIn.goodpix])^2 *  n_elements(*SignalIn.err) / n_elements(err)$
    else $
       minerr = min(*SignalIn.err)^2 *  n_elements(*SignalIn.err) / n_elements(err)
    negative = where(err le minerr, cnt)  
    if cnt ne 0 then err[negative] = minerr  
    if keyword_set(flux_conserve) then *SignalOut.err = sqrt(err) $
    else *SignalOut.err = sqrt(err) / flat 

; When the input pixels are uncorrelated (consecutive pixels independent) and
; we bin to coarser bins (the output will also be independent) we rebin the 
; variance. 
; When we rebin to finer bins, we have to rebin the error (so have to
; divide by sqrt(dof_factor). When the input spectrum is dependent, and we
; go to a coraser grid ... it depends also of the final dof_factor.
    if dof_factor gt 1 then begin  ; the resampling make dependent pixels
       *SignalOut.err /= sqrt(dof_factor)
    endif else if SignalIn.dof_factor gt 1 then begin ; the input is dependent
       d = dof_factor
       if d*SignalIn.dof_factor lt 1 then d=1./SignalIn.dof_factor
       *SignalOut.err /= sqrt(d)
    endif

endif else if n_elements(*SignalIn.err) gt 0 then $
  message, 'Error spectrum and data do not have the same number of pixels.'


if n_elements(*SignalIn.goodpix) ne 0 then begin
; need to (1) make an array for the pixel list, (2) rebin it (3) remake a list
   maskI = bytarr(n_data)
   maskI[*SignalIn.goodpix] = 1B
   maskI = reform(maskI, dim, /OVER)
   maskO = xrebin(borders, maskI, NewBorders, SPLINF=splinf) / flat
   goodpix = where(abs(maskO-1) lt 0.1, cnt)
   if cnt gt 0 then *SignalOut.goodpix = goodpix $
   else message, 'No good pixel left in the output spectrum'
endif

SignalOut.title = SignalIn.title
if n_elements(*SignalIn.hdr) ne 0 then *SignalOut.hdr = *SignalIn.hdr
SignalOut.start = logStart
SignalOut.step =  logScale
SignalOut.sampling = 1
SignalOut.dof_factor = max([1d, SignalIn.dof_factor * dof_factor])

if statmath eq 0 then statmath = check_math(MASK=32)  ; ignore any underflow happening in this routine

sxaddpar, *SignalOut.hdr, 'HISTORY', 'uly_spect_logrebin, vsc='+strtrim(vsc,2)

return, SignalOut

; test the difference between the orig spectrum and the log and
; lin-rebinned one
; The experiment below showed that SPLINF is a little better (and 3-4
; times faster) than SPLINE
;bb =xrebin(borders,*SignalIn.data,NewBorders,/SPLINF)
;ba = xrebin(NewBorders,bb,borders,/SPLINF)
;
;aa = xrebin(NewBorders,*SignalOut.data,borders,/SPLINE)
;
;;print, stddev(*signalIn.data/aa), stddev(*signalIn.data/ba), mean(*signalIn.data/aa), mean(*signalIn.data,/D)/mean(aa,/D)
;print, stddev(*signalIn.data/aa), stddev(*signalIn.data/ba), stddev(*signalIn.data/aa)-stddev(*signalIn.data/ba)

end
