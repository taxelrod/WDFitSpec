;+
; NAME:           
;                 ULY_SPECT_ALLOC, ULY_SPECT_FREE
;
; PURPOSE:       
;                 Allocate/free a spectrum structure
;
; USAGE:
;                 Spec = uly_spect_alloc[, TITLE=title][, DATA=data][, HEADER=header] 
;                                       [, START=start][, STEP=step][, SAMPLING=sampling] 
;                                       [, WAVELEN=wavelen][, ERR=err][, GOODPIX=goodpix] 
;                                       [, DOF_FACTOR=dof_factor]
;                                       [, SPECTRUM=spectrum]
;
;                 uly_spect_free, Spec
;
; KEYWORDS:
;   TITLE:        String describing the spectrum data
;
;   DATA:         Data array
;
;   HEADER:       Header (FITS structure)
;
;   START:        Wavelength of the first pixel (relevant if sampling=0/1)
;
;   STEP:         Step between 2 pixels (relevant if sampling=0/1)
;
;   SAMPLING:     Sampling mode (0: linear, 1: logarithmic, 2: wavelengths list)
;
;   WAVELEN:      Wavelength array, relevant if sampling=2
;
;   ERR:          Error spectrum
;
;   GOODPIX:      Good pixels list
;
;   DOF_FACTOR:   Degree of freedom factor. Ratio between the actual
;                 number of pixels and the number of independent measurements
;                 (i.e. pixels on the detector). This parameter increases
;                 when the spectrum is rebinned to smaller pixels.
;
;   SPECTRUM:     Spectrum structure, initialize the new structure with
;                 the content of this input variable.
;
; DESCRIPTION:
;     Allocate/free a structure for all the spectra with linear or log 
;     (i.e. linear in velocity) wavelength scaling or with explicit wavelength
;     array.
;
;     The different optional keywords that may be passed to this routine
;     indicate how to initialize the structure.
;
;     ULY_SPECT_ALLOC defines an anonymous structure containing the
;     following tags:
;     .title      string
;                 Title                    
;     .hdr        Array of character strings
;                 FITS style header
;     .data       integer, real or double precision array
;                 Pointer on Array of the values of the flux  
;     .err        integer, real or double precision array
;                 Error spectrum
;     .wavelen    Array of double
;                 Used only if 'sampling'=2
;     .goodpix    integer array
;                 List of valid pixels
;     .start      double precision
;                 Wavelength of the centre of the first pixel
;                 If 'sampling'=0, wavelength in Angstrom
;                 If 'sampling'=1, log(wavelength)  (natural log)
;     .step       double precision
;                 Wavelength step between two pixels 
;                 If 'sampling'=0, in Angstrom
;                 If 'sampling'=1, logarithmic step
;     .sampling   integer
;                 Type of sampling scale
;                 0 for linear in wavelength
;                 1 for logarithmic scale  
;                 2 for uneven sampling. Wavelength listed in 'wavelen'
;     .dof_factor Rebinning factor
;
; SEE ALSO:    ULY_SPECT_GET
;
; AUTHOR:      Martin France, Mina Koleva, Philippe Prugniel

;-
; CATEGORY:    ULY
;------------------------------------------------------------------------------
pro uly_spect_free, spectrum

compile_opt idl2
on_error, 0

ptr_free, spectrum.data
ptr_free, spectrum.err
ptr_free, spectrum.wavelen
ptr_free, spectrum.hdr
ptr_free, spectrum.goodpix

end

;------------------------------------------------------------------------------
function uly_spect_alloc, TITLE=title, DATA=data,                    $
                          START=start, STEP=step, SAMPLING=sampling, $
                          ERR=err, WAVELEN=wavelen, GOODPIX=goodpix, $
                          HEADER=header, DOF_FACTOR=dof_factor,      $
                          SPECTRUM=spectIn


compile_opt idl2
on_error, 0


spectrum = { $
             title:'', $
             hdr:ptr_new(/ALLOCATE_HEAP), $
             data:ptr_new(/ALLOCATE_HEAP), $
             err:ptr_new(/ALLOCATE_HEAP), $
             wavelen:ptr_new(/ALLOCATE_HEAP), $
             goodpix:ptr_new(/ALLOCATE_HEAP), $
             start:1.0d, $
             step:1.0d, $
             sampling:1, $
             dof_factor:1. $
           }

if uly_spect_get(spectIn, /VALID) then begin
    if n_elements(*spectIn.hdr) ne 0     then *spectrum.hdr     = *spectIn.hdr
    if n_elements(*spectIn.data) ne 0    then *spectrum.data    = *spectIn.data
    if n_elements(*spectIn.err) ne 0     then *spectrum.err     = *spectIn.err
    spectrum.title      = spectIn.title
    spectrum.start      = spectIn.start
    spectrum.step       = spectIn.step
    spectrum.sampling   = spectIn.sampling
    spectrum.dof_factor = spectIn.dof_factor
    if tag_exist(spectIn, 'wavelen', /TOP_LEVEL) then $
      if n_elements(*spectIn.wavelen) ne 0 then *spectrum.wavelen = *spectIn.wavelen
    if tag_exist(spectIn, 'goodpix', /TOP_LEVEL) then $
      if n_elements(*spectIn.goodpix) ne 0 then *spectrum.goodpix = *spectIn.goodpix $
    else if tag_exist(spectIn, 'mask', /TOP_LEVEL) then begin
        *spectrum.goodpix = where(spectIn.mask eq 1)
    endif
    return, spectrum
endif else if size(spectIn,/TYPE) ne 0 then begin
    message, /INFO, '<spectIn> is not a valid spectrum structure'
    return, 0
endif

if n_elements(title) gt 0 then spectrum.title = title
if n_elements(header) gt 0 then *spectrum.hdr = header
if n_elements(data) gt 0 then *spectrum.data = data
if n_elements(err) gt 0 then *spectrum.err = err
if n_elements(wavelen) gt 0 then begin
    if n_elements(sampling) gt 0 then if sampling ne 2 then $
      message, 'inconsistency, if WAVELEN is passed, SAMPLING must be 2 or ommitted'
    sampling = 2
    *spectrum.wavelen = wavelen
endif
if n_elements(goodpix) gt 0 then *spectrum.goodpix = goodpix
if n_elements(start) gt 0 then spectrum.start = start
if n_elements(step) gt 0 then spectrum.step = step
if n_elements(sampling) gt 0 then spectrum.sampling = sampling
if n_elements(dof_factor) gt 0 then spectrum.dof_factor = dof_factor

return, spectrum

end

;****************************************************************************
