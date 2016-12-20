;+
; NAME:    	ULY_SPECT_GET
;
; PURPOSE:   	Retrieve information from a spectrum structure.
;
; USAGE:
;            	info = uly_spect_get(SignalIn
;                                       [,/WAVERANGE]   
;                                       [,/GOODPIX] 
;                                       [,/HDR]
;                                       [,/MASK]
;                                       [,/VALID], 
;                                       [STATUS=status])
;
; INPUT:
;   SignalIn:
;		Input data structure (spectrum)
;
; KEYWORD:
;     Only one keyword is processed, according to the following priority.
;     If no keyword is provided, the task assumes /VALID .
;
;   [/WAVERANGE]:
;		Returns a 2 elements double precision array giving the
;		wavelength range.
;
;   [/GOODPIX]:
;        	Returns the goodpixels from the signal structre, if not set
;		the task assumes all pixels are good
;     
;   [/HDR]:
;        	Returns the FITS header of the input file, if it exists. If
;		not it returns an empty string. 
;
;   [/MASK]:
;        	Returns a byte array containing the pixel mask (1 for good
;		pixels, 0 for bad).
;		The mask has the same length as the spectrum vector.
;               If the input spectrum structure does not contains a good pixels list, 
;               it is assumed that all the pixels are good.
;
;   [/VALID]:
;        	Returns a boolean: 
;			0 if SignalIn is not a valid spectrum structure;
;                       1 if SignalIn it is a valid structure
;
;   [STATUS=status]
;        	Returnes status code: 0 in case of success, 1 otherwise.
;
; RETURN VALUE:	
;   INFO:
;		See the description of each keyword for a description of
;		the return values. In the case of an error, the return
;		value is -1.
;
; DESCRIPTION:
;     		The input data must be a 'spectrum structure', see
;		ULY_SPECT_ALLOC for more information. 
;		This task allows to retrieve various information about the
;		structure.
;
; HISTORY:    
;		Ph. Prugniel, 2008/01/31, created
;-
; CATEGORY:    ULY
;------------------------------------------------------------------------------

function uly_spect_get, SignalIn, $
                        WAVERANGE=waverange, $
                        GOODPIX=goodpix, $
                        HDR=hdr, $
                        VALID=valid, $
                        MASK=mask, $
                        STATUS=status

if n_params() eq 0 then begin
    message, 'Usage: uly_spect_get, SignalIn, <optional kw>', /INFO
    status = 1
    return, -1
endif

status = 0  ; initalize on 'no error'

valid = 1b
if size(SignalIn, /TYPE) ne 8 then valid = 0b
if valid eq 1 then if ptr_valid(SignalIn.data) eq 0 then valid = 0b
if valid eq 1 then if n_elements(*SignalIn.data) eq 0 then valid = 0b

if keyword_set(waverange) then begin
    if not valid then begin
        status = 1
        return, -1
    endif
    if SignalIn.sampling eq 0 then return, SignalIn.start + $
                [0d,((size(*SignalIn.data))[1]-1)*SignalIn.step]
    if SignalIn.sampling eq 1 then return, exp(SignalIn.start + $
                [0d,((size(*SignalIn.data))[1]-1)*SignalIn.step])
    if SignalIn.sampling eq 2 then begin
        if n_elements(*SignalIn.data) ne n_elements(*SignalIn.wavelen) then $
          return, -1
        return, [(*SignalIn.wavelen)[0], $
                 (*SignalIn.wavelen)[n_elements(*SignalIn.wavelen)-1]]
    endif
endif

if keyword_set(goodpix) then begin
    if not valid then begin
        status = 1
        return, -1
    endif
    if n_elements(*SignalIn.goodpix) ne 0 then $
      goodpix = *SignalIn.goodpix else begin
        npix = (size(*SignalIn.data, /DIM))[0]
        goodpix = lindgen(npix)
    endelse
    return, goodpix
endif

if keyword_set(hdr) then begin
    if not valid then begin
        status = 1
        return, -1
    endif
    hdr = ''
    if tag_exist(SignalIn, 'hdr') then begin
        if ptr_valid(SignalIn.hdr) then $
          if n_elements(*SignalIn.hdr) gt 0 then hdr = *SignalIn.hdr 
    endif 
    return, hdr
endif

if keyword_set(mask) then begin
   if not valid then begin
      status = 1
      return, -1
   endif
   msk = bytarr(n_elements(*SignalIn.data))
   if n_elements(*SignalIn.goodpix) gt 0 then msk[*SignalIn.goodpix] = 1 $
   else msk += 1B
   msk = reform(msk, size(*SignalIn.data, /DIM))
   return, msk
endif

return, valid

end
