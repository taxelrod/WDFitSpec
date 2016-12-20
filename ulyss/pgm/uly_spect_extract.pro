;+
; NAME:      ULY_SPECT_EXTRACT
;
; PURPOSE:   Extract part of a spectrum
;
; USAGE:
;            SignalExtr = uly_spect_extract(SignalIn[, ONED=pos] 
;                        [, WAVERANGE=waverange][,/OVERWRITE][,STATUS=status])
; ARGUMENTS:
;     SignalIn  input data structure (spectrum)
;
; INPUT KEYWORD:
;     ONED   Position of a 1D scan to extract from the file
;            Scalar or array giving the position of the scan to extract
;            along each axis after the first.
;
;     WAVERANGE A 1 or 2 element array specifying the wavelength region to 
;            extract. It must be ordered [low,high]. If only one element
;            is given, it specifies the minimum wavelength to extract.
;
;     OVERWRITE Set this keyword if SignalExtr is the same variable as
;            SignalIn to prevent a memory leak
;
; RETURN VALUE:
;     1D extracted spectrum (spect structure, see uly_spect_alloc)
;
; OUTPUT KEYWORD:
;     STATUS=status
;        status code: non-null in case of error
;
; DESCRIPTION:
;     Extract part of a spectrum and return a 'spect' structure.
;
;     Two options may be specified separately or together:
;     - ONED is used to extract a 1D spectrum from a 2D dataset,
;            like for example a long-slit spectrum.
;            If the specified position is invalid, an informational
;            message is issued, status is set to 1, and the routine
;            returns -1.
;     - WAVERANGE is used to extract a wavelength region. 
;            The output spectrum is cut to the largest range containing
;            the specified limits.
;            If the specified value is not a 1 or 2 elements array, or if
;            the specification has no intersection the the actual
;            range, and informational message is issued, status is set
;            to 1, and -1 is returned.
;
;     The input data must be a 'spectrum structure', see ULY_SPECT_ALLOC.
;
; HISTORY:    Ph. Prugniel, 2008/01/31, created
;-
; CATEGORY:    ULY
;------------------------------------------------------------------------------
function uly_spect_extract, SignalIn, ONED=pos, WAVERANGE=waverange, $
                            OVERWRITE=overwrite, STATUS=status

  if n_params() ne 1 then begin
     message, 'Usage: s=uly_spect_extract(SignalIn, ONED=pos, WAVERANGE=waverange)', /INFO
     status = 1
     return, -1
  endif

  status = 0                    ; initialize on 'no error'

  if not uly_spect_get(SignalIn, /VALID) then begin
     message, '<SignalIn> is not a valid spectrum', /INFO
     status = 1
     return, -1
  endif

  if keyword_set(overwrite) then SignalOut=SignalIn else $
     SignalOut = uly_spect_alloc(SPECTRUM=SignalIn)

  
; extraction of a 1D scan
  ndim = size(*SignalOut.data, /N_DIM)
  if n_elements(pos) ne 0 and ndim eq 1 then begin
     if pos[0] ne 0 then begin
        message, 'Invalid ONED position', /INFO
        uly_spect_free, SignalOut
        status = 1
        return, -1
     endif
  endif else if n_elements(pos) ne 0 then begin
     ndim = size(*SignalOut.data, /N_DIM)
     if n_elements(pos) ge ndim and ndim gt 1 then begin
        message, 'Invalid ONED position', /INFO
        uly_spect_free, SignalOut
        status = 1
        return, -1
     endif

     dim = size(*SignalOut.data, /DIM)
     n = where(pos lt 0 or pos ge dim[1:*], cnt)
     if cnt gt 0 then begin
        message, 'Invalid ONED position', /INFO
        uly_spect_free, SignalOut
        status = 1
        return, -1
     endif

     sdim = [1, exp(total(alog(dim[1:*]), /CUMUL))]
     ntot = n_elements(*SignalOut.data)
     ind = total([pos]*sdim)

     *SignalOut.data = (reform(*SignalOut.data, dim[0], ntot/dim[0]))[*,ind]
     if n_elements(*SignalOut.err) eq ntot then $
        *SignalOut.err = (reform(*SignalOut.err, dim[0], ntot/dim[0]))[*,ind]

     msk = uly_spect_get(SignalOut, /MASK)
     if n_elements(msk) eq ntot then begin
        msk = (reform(msk,ntot/dim[0]))[*,ind]
        *SignalOut.goodpix = where(msk eq 1)
     endif

  endif

; cut a wavelength region
  if n_elements(waverange) ne 0 then begin
     if n_elements(waverange) gt 2 then begin
        message, 'WAVERANGE must be a 1 or 2 elements array', /INFO
        uly_spect_free, SignalOut
        status = 1
        return, -1
     endif
     range = uly_spect_get(SignalOut,/WAVERANGE) 
     if range[1] lt waverange[0] then begin
        message, 'WAVERANGE ('+string(waverange,FOR='(2f)')+ $
                 ') is out of bounds ('+string(range,FOR='(2f,1H))'), /INFO
        uly_spect_free, SignalOut
        status = 1
        return, -1
     endif
     if n_elements(waverange) eq 2 then $
        if waverange[0] gt waverange[1] then begin
        message, 'WAVERANGE must be ordered low<high', /INFO
        uly_spect_free, SignalOut
        status = 1
        return, -1
     endif
     if n_elements(waverange) eq 2 then $
        if range[0] gt waverange[1] then begin
        message, 'WAVERANGE is out of bounds', /INFO
        uly_spect_free, SignalOut
        status = 1
        return, -1
     endif
     
     if SignalOut.sampling eq 0 or SignalOut.sampling eq 1 then begin
        wr = waverange
        if SignalOut.sampling eq 1 then wr = alog(waverange)
        
        nummin = floor((wr[0]-SignalOut.start) / SignalOut.step + 0.01d)
        if nummin lt 0 then nummin = 0
        
        npix = (size(*SignalOut.data, /DIM))[0]
        if n_elements(waverange) eq 2 then begin
           nummax = ceil((wr[1]-SignalOut.start) / SignalOut.step - 0.01d)
           if nummax ge npix then nummax = npix-1
        endif else nummax = npix-1
        
        SignalOut.start += nummin * SignalOut.step

        if n_elements(*SignalOut.goodpix) ne 0 then begin 
           msk = (uly_spect_get(SignalOut, /MASK))[nummin:nummax, *]
           *SignalOut.goodpix = where(msk eq 1)
        endif
        s = size(*SignalOut.data, /DIM) ; data and err may be more than 1D
        *SignalOut.data = reform(*SignalOut.data, s[0], n_elements(*SignalOut.data)/s[0], /OVER)
        s[0] = nummax - nummin + 1
        *SignalOut.data = reform((*SignalOut.data)[nummin:nummax,*], s, /OVER)
        
        if n_elements(*SignalOut.err) ne 0 then begin
           s = size(*SignalOut.err, /DIM)
           *SignalOut.err = reform(*SignalOut.err, s[0], n_elements(*SignalOut.err)/s[0], /OVER)
           s[0] = nummax - nummin + 1
           *SignalOut.err = reform((*SignalOut.err)[nummin:nummax,*], s)
        endif
     endif else if SignalOut.sampling eq 2 then begin
        if n_elements(waverange) eq 0 then lmn = min(*SignalOut.wavelen) $
        else lmn = waverange[0]
        if n_elements(waverange) le 1 then lmx = max(*SignalOut.wavelen) $
        else lmx = waverange[1]

        extr = where(*SignalOut.wavelen ge lmn and *SignalOut.wavelen le lmx, cnt)
        dim = size(*SignalOut.data, /DIM) ; data and err may be more than 1D
        nscan =  n_elements(*SignalOut.data)/dim[0] ; number of 1D spectra
        
        if cnt eq 0 then begin
           message, 'No data left in extraction', /INFO
           uly_spect_free, SignalOut
           status = 1
           return, -1
        endif

        *SignalOut.wavelen = (*SignalOut.wavelen)[extr,*]
        if n_elements(*SignalOut.goodpix) gt 0 then begin
            mask = bytarr(n_elements(*SignalOut.data))
            mask[*SignalOut.goodpix] = 1
            mask = reform(mask, dim[0], nscan)
            mask = mask[extr,*]
            *SignalOut.goodpix = where(mask eq 1)
        endif 
        *SignalOut.data = reform(*SignalOut.data, dim[0], nscan, /OVER)
        *SignalOut.data = (*SignalOut.data)[extr,*]
        if n_elements(*SignalOut.err) ne 0 then begin
            *SignalOut.err = reform(*SignalOut.err, dim[0], nscan, /OVER)
            *SignalOut.err = (*SignalOut.err)[extr,*]
        endif
        dim[0] = (size(*SignalOut.data, /DIM))[0]
        *SignalOut.data = reform(*SignalOut.data, dim, /OVER)
        if n_elements(*SignalOut.err) ne 0 then $
           *SignalOut.err = reform(*SignalOut.err, dim, /OVER)
        
     endif
     
  endif

  return, SignalOut

end
