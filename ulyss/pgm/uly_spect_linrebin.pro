;+
; NAME:      
;                 ULY_SPECT_LINREBIN
;
; PURPOSE:   
;                 Rebin a spectrum in linear wavelength scale
;
; USAGE:     
;                 spect = uly_spect_linrebin(<SignalIn>, <step>)
;
; ARGUMENTS:
;   <SignalIn>:   Input data structure (spectrum)
;
;   <Step>:       Size of each pixel in Angstrom for the output data
;
; KEYWORDS:
;   WAVERANGE     Input scalar or 2-elements vector. Wavelength range.
;                 If a single value s passed, it is interpreted as the
;                 start of the wavelength range.
;                 By default, use the current range.
;
;   FLUX_CONSERVE Input boolean.
;                 Set this keyword to skip the intensity correction
;                 for the change of the local pixel size.
;                 If flux is conserved, intensity is not; by default
;                 intensity is conserved.
;
;   EXACT         Input boolean.
;                 Specify that the output must match exactly the WCS
;                 specified by <waverange> and <vsc>. Otherwise approximations
;                 can be made to avoid shifts by small fractions of pixels.
;
;   OVERWRITE     Input boolean.
;                 Set this keyword if input and output spect are the
;                 same variable as to prevent a memory leak
;
; OUTPUT:
;                 Spectrum structure containing the rebinned data
;
; DESCRIPTION:
;     This routine performs rigorously a rebinning into linearly
;     spaced wavelength points. The algorithm is similar to ULY_STECT_LOGREBIN.
;
;     The input data must be a 'spectrum structure' (see ULY_SPECT_ALLOC), 
;     and the data can be sampled linearly, in log or at explicit wavelength.
;
;     If <step> is initialized in input, this value will be used to rebin
;     the input spectrum; otherwise this step is computed in order to keep 
;     the same number of pixels in output as in input,the computed
;     value is returned in <step>.
;
;     The error spectrum and pixel mask, if present in (SignalIn) are
;     also transformed in accordance with the signal. 
;
;     If the data and error arrays contained in (SignalIn) have two dimensions,
;     each line is rebinned similarly. 
;
; REQUIRED ROUTINES:  
;                 XREBIN.PRO
;
; HISTORY:        Mina Koleva, 2007/03/01, created

;-
; CATEGORY:       ULY
;------------------------------------------------------------------------------
function uly_spect_linrebin, SignalIn, step, WAVERANGE=waverange, $
                             FLUX_CONSERVE=flux_conserve,         $
                             EXACT=exact, OVERWRITE=overwrite

compile_opt idl2
on_error, 2

if size(SignalIn,/TYPE) ne 8 then message, 'Input spectrum is undefined'

if n_elements(*SignalIn.data) eq 0 then message, 'Input spectrum undefined'

s = size(*SignalIn.data)
npix = (size(*SignalIn.data,/DIM))

if SignalIn.sampling eq 0 then begin
    cpout = 0
;    Note that comparison between SignalIn.step and step should tolerate
;    rounding errors, this is not the case yet.
    if (n_elements(step) eq 1) then begin
        if SignalIn.step eq step then cpout = 1
    endif else if n_elements(step) eq 0 then cpout = 1
    if keyword_set(exact) and cpout eq 1 and n_elements(waverange) gt 0 then begin
        nshift = (SignalIn.start-waverange[0])/SignalIn.step mod 1d
        if abs(nshift) gt 0.001 then cpout = 0

        if cpout eq 1 and n_elements(waverange) eq 2 then begin
           nshift = (SignalIn.start + (npix[0]-1)*SignalIn.step $
                    - waverange[0])/SignalIn.step mod 1d
           if abs(nshift) gt 0.001 then cpout = 0
        endif
    endif
    if cpout eq 1 then begin    ; just have to copy input to output...
        step = SignalIn.step
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


; Determine the step, preserving the total number of pixels
;  (this block is useless as long as we do not have the WAVERANGE kw ..)
if n_elements(step) eq 0 then begin
    if SignalIn.sampling eq 0 then step = SignalIn.step $
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
            nw = (alog(waverange) - SignalIn.start)/SignalIn.step -0.5d
            if n_elements(waverange) gt 1 then $
              nw[1] = min([nw[1]+1, npix[0]-0.5]) $
            else nw = [nw[0], npix[0]-0.5]
        endif else nw = [-0.5d, npix[0]-0.5]

        wr = exp(SignalIn.start + nw * SignalIn.step)
        step = (wr[1]-wr[0]) / (nw[1] - nw[0])
    endelse
endif

; Determine the start wavelength, linStart, and number of pix, nout, in output
linRange = SignalIn.start + [-0.5d,(npix[0]-0.5)] * SignalIn.step
if SignalIn.sampling eq 1 then linRange = exp(linRange) else $
if SignalIn.sampling eq 2 then linRange = [(*SignalIn.wavelen)[0],(*SignalIn.wavelen)[npix-1]]
linRange += [0.5d,-0.5d]*step
if n_elements(waverange) gt 0 then begin
   nshift = ceil(max([0d, (linRange[0] - waverange[0]) / step]))
   linRange[0] = waverange[0] + step * nshift  ; center of 1st out pix

   if n_elements(waverange) eq 2 then $
      linRange[1] = min([waverange[1], linRange[1]]) 

   if linRange[1] lt linRange[0] then begin
       message, /INFO, 'waverange is not valid: '+ $
         strjoin(strtrim(waverange,2),',')
       return, 0
   endif
endif 
nout = round((linRange[1]-linRange[0]) / step + 1)
linStart = linRange[0]

; determine the DOF factor
if SignalIn.sampling lt 2 then begin
   if SignalIn.sampling eq 1 then linRange = alog(linRange)
   nin = round((linRange[1]-linRange[0]) / SignalIn.step + 1)
endif else nin = n_elements(*SignalIn.wavelen)
dof_factor = nout/double(nin)

if linStart-step/2d gt borders[npix[0]] then begin
    message, 'start value is not in the valid range', /INFO
    return, 0
endif

; define the new borders of the bins
NewBorders = linStart + (dindgen(nout+1)-0.5d) * step

; the rest of this routine is exactly the same as uly_spect_logrebin
; (shall we put it in a common low-level routine?)

dim = size(*SignalIn.data, /DIM)
n_data = n_elements(*SignalIn.data)
n_err = n_elements(*SignalIn.err)
n_dim = size(*SignalIn.data, /N_DIM)

; Determine the conversion factor/vector accounting for pixel size change
case SignalIn.sampling of
    0: flat = step/SignalIn.step
    1: begin        
        flat = step / SignalIn.step / (linStart + dindgen(nout)*step)
        if n_dim gt 1 then flat = rebin(flat, [nout, dim[1:*]])
    end
    else: begin
        flat = fltarr(dim[0], /NOZERO) ; faster to create an array and
        replicate_inplace, flat, 1B    ; populate with 1s 
        flat = xrebin(borders, flat, NewBorders, /SPLINF) 
        if n_dim gt 1 then flat = rebin(flat, [nout, dim[1:*]])
    end
endcase


if keyword_set(overwrite) then SignalOut = SignalIn $
else SignalOut = uly_spect_alloc()

SignalIn_type = size(*SignalIn.data, /TYPE)

if keyword_set(flux_conserve) then $
  *SignalOut.data = xrebin(borders,*SignalIn.data,NewBorders,/SPLINF) $
else $
  *SignalOut.data = xrebin(borders,*SignalIn.data,NewBorders,/SPLINF) / flat

if SignalIn_type lt 5 then *SignalOut.data = float(temporary(*SignalOut.data))

if n_err eq n_data then begin
; err is very small, then the ^2 is making it close to 0 and log_rebinning
; may produce negative values.
; This cause errors "% Program caused arithmetic error: Floating
; illegal operand", and leave NaN values in the err spectrum. 
; To prevent this we set these negative values to the min of the positive.
    err = xrebin(borders,*SignalIn.err^2,NewBorders,/SPLINF) 
    if n_elements(*SignalIn.goodpix) ne 0 then $
       minerr = min((*SignalIn.err)[*SignalIn.goodpix])^2 *  n_elements(*SignalIn.err) / n_elements(err)$
    else $
       minerr = min(*SignalIn.err)^2 *  n_elements(*SignalIn.err) / n_elements(err)
    negative = where(err le minerr, cnt)  
    if cnt ne 0 then err[negative] = minerr 
    if keyword_set(flux_conserve) then *SignalOut.err = sqrt(err) $
    else *SignalOut.err = sqrt(err) / flat

; case of correlated pixels, see uly_spect_logrebin
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
; need to (1) make a spectrum for the pixel list, (2) rebin it (3)
; remake a list
    maskI = replicate(0, npix[0])
    maskI[*SignalIn.goodpix] = 1
    maskO = xrebin(borders, maskI, NewBorders, /SPLINF) / flat
    *SignalOut.goodpix = where(abs(maskO-1) lt 0.1)
endif

SignalOut.title = SignalIn.title
if n_elements(*SignalIn.hdr) ne 0 then *SignalOut.hdr = *SignalIn.hdr
SignalOut.start = linStart
SignalOut.step =  step
SignalOut.sampling = 0  ; linear sampling
SignalOut.dof_factor = max([1d, SignalIn.dof_factor * dof_factor])

sxaddpar, *SignalOut.hdr, 'HISTORY', 'uly_spect_linrebin, step='+strtrim(step,2)
return, SignalOut

end

pro uly_test_spect_rebin

; a spectrum used for test
common uly_path, uly_root
if n_elements(uly_root) ne 0 then begin
    if file_test(uly_root+'/data/VazMiles_z-0.40t07.94.fits') eq 1 then $
      uly_datadir = uly_root+'/data'
endif

if n_elements(uly_datadir) ne 0 then $
  filetest = uly_datadir + '/VazMiles_z-0.40t07.94.fits' $
else $
  message, 'Could not find the test file'

lmin =  [4000.,5905.,6300.]
lmax =  [5880.,6260.,6800.]

; read the spectrum (the output is log_rebinned...
spectrum = uly_spect_read(filetest, lmin, lmax)
spec_log = uly_spect_logrebin(spectrum, vsc2)
spec_lin = uly_spect_linrebin(spec_log, step)

relat_dif = (*spectrum.data - *spec_lin.data) / *spectrum.data

print, 'n_elem (in) =', n_elements(*spectrum.data), $
  ' (out) =', n_elements(*spec_lin.data)

print, 'mean (in) =', mean(*spectrum.data), ' (out) =', mean(*spec_lin.data)

print, 'relat resid = ', mean(relat_dif)

if abs(mean(relat_dif)) lt 1d-3 then print, 'Test is sucessful' $
else print, 'Test failed'

end
