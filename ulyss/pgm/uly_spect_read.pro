;+
; NAME:
;              ULY_SPECT_READ
;
; PURPOSE:
;              Read a a spectrum from a FITS file
;
; USAGE:
;    spec = uly_spect_read(<file_in>, 
;                          <lmin>, <lmax>
;                          [, VELSCALE=<velscale>]
;                          [, SG=<sg>]
;                          [, ERR_SP=<err_sp>]
;                          [, SNR_SP=<snr_sp>]
;                          [, MSK_SP=<msk_sp>]
;                          [, DISP_AXIS=<disp_axis>]
;                          [, /QUIET] )
;
; ARGUMENTS:
;    file_in: 
;	Name or logical unit (lun) of the FITS file containing a spectrum.
;	Passing the lun of an already open file is particularly useful to
;	read data from a given FITS extension (i. e. not necessarly the
;	first). 
;
;       For example:
;          lun = fxposit('file.fits', 1)  ; to skip the first extension
;          spec = uly_spect_read(lun) ; to read from the 2nd extension
;
;    lmin,lmax: (input)
;       Wavelength range to use (maybe several intervals)
;                                                            
; KEYWORDS:
;    VELSCALE:  Input scalar. Optional.
;               Velocity scale (km/s), i.e. size of the pixel after 
;               log-rebinning.
;               If specified call ULY_SPECT_LOGREBIN to resample the spectrum.
;               determined with ULY_SPECT_LOGREBIN.
;
;    ERR_SP:    Input character string or integer. Optional.
;               Name or lun of the input file containing the error spectrum.
;               This is used in case the error spectrum is not contained
;               in <file_in>. The error spectrum must have the same
;               WCS as <file_in>.
;
;    SNR_SP:    Input character string or integer. Optional. 
;               Name or lun of the input file containing the SNR spectrum.
;               Used if the error spectrum is not contained in <file_in>
;               and if ERR_SP is not specified.
;
;    MSK_SP:    Input character string or integer. Optional. 
;               Name or lun of the input file containing the mask spectrum.
;               Used if the mask is not contained in <file_in>. The mask
;               must have the same WCS as <file_in>.
;
;    SG:        Input. 
;               Redshift (dimensionless) used to change the spectral WCS from
;               the current frame (probably heliocentric) to the rest-frame.
;
;    FORMAT:    Input integer (1, 2 or 3). Optional    
;               Bypass the automatic format recognition and force the format.
;               See below for description of the formats.
;
;    DISP_AXIS: Input. 
;               Number of the wavelength axis. Default = 1; Valid values 1/2.
;               This information is used only for the format 3 (see the
;               section 'Description'), and superseeds indications that 
;               would be read from the header of the FITS file.
;
;    QUIET:     Input
;               Verbosity control. By default the function prints
;               some information about the data it reads.
;
; OUTPUT:
;    spect structure, see ULY_SPECT_ALLOC
;
; DESCRIPTION:
;   Read a spectrum from a FITS file, including possibly the associated
;   noise and a mask to exclude the bad regions, shift it back to the 
;   restframe (if SG is given) and rebin it to a logarithmic
;   wavelength scale (if VELSCALE is given).
;
;   The return value, <spec>, is a structure containing pointers
;   and has to be freed in the calling procedure as: uly_spect_free, spec
;   See the documentation of ULY_SPECT_ALLOC for a description of
;   this structure.
;
;   Three different formats of FITS files for spectra are supported:
;   1- SDSS style 
;   2- Binary tables (BINTABLE)
;   3- Image type array for 1D and long-slit spectra
;
;   In all cases, the spectrum is supposed to contain fluxes, not log(flux).
;   Whether they are absolutely, relatively or not flux-calibrated, or 
;   normalized to a pseudo-continuum, is not coded in the output structure
;   and is usually not critical for the analysis performed by the package. 
;   Note however that the original FITS header is attached in the output
;   structure, so, if necessary this information can be retrieved after
;   reading the data.
;
;   Description of the SDSS style (< DR9), format 1:
;     The SDSS style allows to store wavelength resampled spectra
;     together with their associated error, mask of bad regions, and
;     potentially any other additional data, like for example a best
;     fitted model (see our format for the solution of the fit, 
;     described in the documentation of ULY_SOLUT_READ).
;
;     The various quantities, spectrum and associated data, are stored
;     on the last axis of an image array and described by the keywords
;     ARRAYi, where i is the index along the last axis. The following
;     values of are used by ULY_SPECT_READ:
;        DATA      Flux spectrum (SPECTRUM is an alternative)
;        ERROR     Error spectrum (1 sigma)
;        MASK      Mask of good pixels (0=bad)
;
;     The ARRAYi keywords are used to recognize the SDSS format. Use the FORMAT
;     keyword (see above) if those keywords are absent. 
;
;     Although this format is used by the SDSS project for 1D spectra only,
;     it is easily extended here for larger dimensions by assuming that
;     the associated data (DATA, ERROR...) are stacked on the last axis.
;
;     Note also that in the SDSS project, the spectral WCS is not described
;     by a standard WCS and is difficult to automatically guess from
;     the other keywords (inheritance of keywords along the various steps
;     of the processing gives room to ambiguous interpretation).
;     The log10 scaling of SDSS data is recognized by the automatic process
;     described below (though standard wcs should be in natural log).
;
;     Note that since DR6 in SDSS the "DATA" array is called "SPECTRUM"
;
;   Description of BINTABLE, format 2:
;     There are many ways to store spectra in FITS binary tables. 
;     The two main options are either to store them as columns of a
;     table, or in a single cell. This latter format, which is not supported
;     here, allows to store the various associated data (arrays or scalars)
;     in a single row of a table, and therefore store several spectra
;     in the file.
;     The present reader handle spectral information stored as columns.
;     Its advantage is to allow to store spectra which are not
;     resampled in wavelength (or log of wavelength): A dedicated column
;     may provide the wavelength of each pixel.
;
;     The current version is limited to the reading of the output of:
;     - the DEIMOS pipeline (spec2d).
;     The following columns are read:
;        SPEC    1D spectrum
;        LAMBDA  wavelength in angstroms.
;        IVAR    inverse variance array
;        ORMASK  mask with 0-good and 1-bad
;     - the SDSS DR9 pipeline.
;        FLUX    1D spectrum
;        LOGLAM  log10 wavelength in angstroms in vacuum, transformed
;        to air and log while reading
;        IVAR    inverse variance array
;        OR_MASK is not read since it does not consist of 0 and 1 
; 
;   Description of image type arrays for 1/2D spectra, format 3:
;      The most basic way to store 1D or long-slit spectra is as a 1/2D array. 
;      The axis (or axes) being described by a WCS.
;
;      The main disavantage of this format is the difficulty to
;      associate the other data (error or mask may be stored in other
;      extensions described in the main header, but there is not
;      standard for this process).
;
;      Use the ERR_SP, SNR_SP and MSK_SP keywords can be used to read the noise
;      and mask for this format.
;
;   Note about the spectral WCS:
;     For the formats 1 and 3, the routine first attempts to read the standard 
;     WCS described by the keywords CTYPEi, CRPIXi, CRVALi and CDi_i or CDELTi, 
;     where i is the number of the spectral axis previously determined.
;
;     If CTYPEi is absent but the other keywords present, the routine supposes
;     that the wavelengths are expressed in Air (not in vacuum) and guess the 
;     wavelength scaling from the world coordinate of the reference pixel 
;     (CRVALi). If this value is less that 4, it is assumed that the axis is 
;     in log10(wave/Angstrom). Between 4 and 9 in ln(wave/Angstrom), like if 
;     CTYPEi='AWAV-LOG', and above that it is the wavelength in Ansgtrom 
;     like if CTYPEi='AWAV'.
;
;     CUNITi is never used (at present).
;
;   Note about the conversion of vacuum wavelengths:
;     The IAU definition for the conversion between Air and Vacuum
;     wavelength, respectively (VAC) and (AIR) is:
;       VAC = AIR * 1 + 6.4328D-5 + 2.94981D-2/(146 - sigma2) + 
;             2.5540D-4/(41 - sigma2)
;       where sigma2 = 1/AIR^2 (the wavelengths are in Angstrom)
;     This formula is cited in Morton 1991 ApJS 77, 119
;
;     So, approximately in the visible range: AIR = VAC / 1.00028
;     (i.e. a shift of 84 km/s).
;     ULY_SPECT_READ applies this approximate conversion to the WCS, to avoid
;     a resampling. The IDL function VACTOAIR shall be used if a higher
;     precision is required.
;
;     The wavelength calibration of the SDSS spectra is in Heliocentric 
;     VACUUM wavelength.
;
; EXAMPLE:
;     Read a spectrum, and plot it:
;       s = uly_spect_read(uly_root+'/data/VazMiles_z-0.40t07.94.fits')
;       uly_spect_plot, s
;     (note that the file name could have been pass directly to uly_spect_plot)
;
; AUTHOR:
;              Mina Koleva, Philippe Prugniel
;
; HISTORY:
;             Yue WU, 2013/01/25 - reading SDSS DR9 format
;             Mina Koleva 2013/03/01 - debuging, cleaning the reading
;                                      of DR9 format                                                           
;-
; CATEGORY:    ULY
;------------------------------------------------------------------------------
function uly_spect_filefmt, hdr, FILE=file
; Analyse the header of a FITS file to guess what format it is
; The routine is presently recognizing 3 spectral formats.
; The aim is perform minimum tests in order to select the proper reading 
; algorithm. It does not warrant that the file is actually readable.
;
; The function returns:
;   0 : The format is not recognized
;   1 : SDSS-like format
;   2 : Spectrum in columns of a BINTABLE
;   3 : 1D or LSS-2D spectrum with a 'fair' WCS
; If the variable hdr does not contain a FITS header, the header can 
; be read from <file>

  if size(hdr,/TYPE) ne 7 or n_elements(hdr) le 1 then begin
     fits_read, file, data, hdr, MESS=mess
     if mess ne '' then begin
        message, mess, /INFO
        return, 0
     endif
  endif

; Search if it is a SDSS-like format (ie. have ARRAYi)
  if sxpar(hdr, 'NAXIS') gt 1 then begin  
     array = sxpar(hdr, 'ARRAY*', COUNT=cnt1)
     naxis = sxpar(hdr, 'NAXIS*', COUNT=cnt2)
     ns = sxpar(hdr, 'NAXIS'+strtrim(cnt2,2))
     if cnt1 gt 0 and ns ge cnt1 then return, 1
  endif

; Search if it is a BINTABLE format 
  if strtrim(sxpar(hdr, 'XTENSION')) eq 'BINTABLE' then begin
     return, 2
  endif

; Search if it can be a 1D or LSS spectrum in image array
  if sxpar(hdr, 'NAXIS') le 3 then begin  
     return, 3
  endif

  return, 0
end

;==============================================================================
; reading routine for SDSS-style format (format=1)
function uly_spect_read_sdss, data, h, ERR_SP=err_sp, SNR_SP=snr_sp, QUIET=quiet

if not keyword_set(quiet) then message, /INFO, 'SDSS style (< DR9)'

; pb with sdss style: the definition of MASK in SDSS original data
; is complex (this is a bit mask, 0 is good)
; We would like to understand also a simpler mask made of 0 and 1s, 1
; is good...

; identify the subarrays that we will extract (the others will be ignored)
array = strtrim(sxpar(h, 'ARRAY*', COUNT=ca), 2)
if ca gt 0 then begin
   n1 = where (array eq 'DATA'  or array eq 'SPECTRUM')
   n2 = where (array eq 'ERROR')
   n3 = where (array eq 'MASK')
   n4 = where (array eq 'WAVELEN')
endif else begin  ; No ARRAY kw, we assume a standard structure
   naxi = sxpar(h, 'NAXIS')
   naxl = sxpar(h, 'NAXIS'+strtrim(naxi,2))
   if naxl gt 0 then n1 = 0 else n1 = -1  ; spectrum array
   if naxl gt 1 then n2 = 0 else n2 = -1  ; error array
   n3 = -1 
   n4 = -1
   if not keyword_set(quiet) then begin
      message, /INFO, 'There is no ARRAYn keywords in the header, assume that the first line is DATA and second ERROR'
   endif
endelse

; analyse the WCS
sampling = -1
vacuum = 0     ; Air or vacuum wavelengths

if n4 lt 0 then begin
;   check CTYPE
    ctype = sxpar(h,'CTYPE1')
    if strmid(ctype,0,4) eq 'WAVE' then vacuum = 1
    if sxpar(h, 'VACUUM') eq 1 then vacuum = 1  ; Used in SDSS

    if strmid(ctype,0,4) eq 'WAVE' or strmid(ctype,0,4) eq 'AWAV' then begin
        if strmid(ctype,5,3) eq 'WAV' or strmid(ctype,5,3) eq '   ' then begin
            sampling = 0
        endif else if strmid(ctype,5,3) eq 'LOG' then begin
            sampling = 1
        endif
    endif

    step = double(string(sxpar(h,'CD1_1',COUNT=count)))
    if count eq 0 then begin
        step = double(string(sxpar(h,'CDELT1',COUNT=count)))
    endif
    if count eq 0 then message, 'Cannot decode the WCS'
    crpix = double(string(sxpar(h,'CRPIX1',COUNT=count))) - 1d
    if count eq 0 then crpix = 0d
    start = double(string(sxpar(h,'CRVAL1', COUNT=count))) - crpix*step
    if count eq 0 then message, 'Cannot decode the WCS'

    if sampling lt 0 then begin
        if start le 4 then begin
            sampling = 1 
            start *= alog(10d)
            step *= alog(10d)
            if not keyword_set(quiet) then $
              message, /INFO, 'Assume that the sampling is in log10'
        endif else if start le 9 then begin
            sampling = 1
            if not keyword_set(quiet) then $
              message, /INFO, 'Assume that the sampling is in ln'
        endif else begin
            sampling = 0
            if not keyword_set(quiet) then $
              message, /INFO, 'Assume that the sampling is linear in wavelength'
        endelse
    endif
    if vacuum eq 1 then begin
       if not keyword_set(quiet) then $
           message, /INFO, 'Wavelength in VACUUM ... approximately converted'
       if sampling eq 0 then begin 
           start /= 1.00028d 
           step /= 1.00028d 
       endif else if sampling eq 1 then start -= 0.00028D
    endif
endif else sampling = 2

; reformat the data array
ndim = size(data,/N_DIM)
dim = size(data, /DIM)

narray = dim[ndim-1]
dim = dim[0:ndim-2]
tot = product(dim)
data = reform(data, tot, narray, /OVER)

spect = uly_spect_alloc(START=start, STEP=step, SAMP=sampling, HEADER=h)

; Remove some WCS keywords  (and ARRAY*) that may be outdated
sxdelpar, *spect.hdr, ['VACUUM', 'CTYPE1', 'CRVAL1', 'CDELT1', 'CD1_1', 'CRPIX1', 'DC-FLAG', 'WAT0_*', 'WAT1_*', 'WFITTYPE', 'ARRAY'+strtrim(indgen(n_elements(array))+1,2), ['NAXIS', 'CTYPE', 'CUNIT'] + strtrim(ndim,2)]

if n1 ge 0 then *spect.data = reform(data[*,n1], dim)
if n2 ge 0 then *spect.err = reform(data[*,n2], dim)

if n3 ge 0 then begin
    m = where(finite(data[*,n3]) eq 1, cnt)
    if cnt gt 0 then begin
        m = where(data[*,n3] ne 0 and data[*,n3] ne 1, cnt)
        if cnt gt 0 then begin
            message, /INFO, $
              'The mask is not made of 0 and 1s ... it is ignored' 
            message, /INFO, '   (a standard mask has 0=bad, 1=good)' 
        endif else $
          *spect.goodpix = where(data[*,n3] gt 0) ;goodpix is a 1D list
    endif
endif

if n4 ge 0 then *spect.wavelen = reform(data[*,n4], dim)

if n_elements(*spect.err) gt 0 then begin
    if n_elements(*spect.goodpix) eq 0 then $
      m = where(*spect.err le 0, cnt, COMP=g) $
    else m = where((*spect.err)[*spect.goodpix] le 0, cnt, COMP=g)
    if cnt gt 0 then begin
        if not keyword_set(quiet) then $
           message, /INFO, 'Pixels with 0 or negative errors were masked'
        if n_elements(*spect.goodpix) eq 0 then *spect.goodpix = g $
        else *spect.goodpix = [*spect.goodpix, g]
    endif
endif

dof_factor = double(string(sxpar(h, 'DOF_FACT', COUNT=count)))
if count eq 1 then spect.dof_factor = dof_factor

return, spect

end

;==============================================================================
; reading routine for BINTABLE format (format=2)
function uly_spect_read_tbl, data, h, ERR_SP=err_sp, SNR_SP=snr_sp, QUIET=quiet
  
; we should there find what columns to read (for data, error and mask)
  
; for the moment it is only suited to the spectra produced by spec2d,
; the pipeline of DEIMOS: 
; http://astro.berkeley.edu/~cooper/deep/spec2d/primer.html

  if tag_exist(data, 'SPEC') then flux = data.spec $
  else if tag_exist(data, 'FLUX') then flux = data.flux 

; in DEIMOS and SDSS DR9: IVAR is the inverse of the variance of each px
  if tag_exist(data, 'IVAR') then begin
     zeros = where(data.ivar eq 0, cnt, compl=goodpix)
     if cnt gt 0 then data.ivar[zeros] = 1 ;we put arbitary value for the masked pxs
     err = 1d/sqrt(data.ivar)
  endif $
  else if tag_exist(data, 'error') then err = data.error 

;; DEIMOS
  if tag_exist(data, 'ORMASK') then begin
     m = where(finite(data.ormask) eq 1, cnt)
     goodpix = where(data.ormask eq 0) 
     message, /INFO, $
              'DEIMOS Mask 0=good, 1=bad' 
  endif
  
  if tag_exist(data, 'LAMBDA') then wave = data.lambda 
  if tag_exist(data, 'WAVE') then wave = data.wave 
  if tag_exist(data, 'LOGLAM') then begin ; SDSS DS9
     if not keyword_set(quiet) then begin
        message, /INFO, 'SDSS DR9: (1) Sampling is in log10, (2) Wavelength in VACUUM ... approximately converted'
     endif
     start = (data.loglam)[0] * alog(10d)  - 0.00028D
     step = ((data.loglam)[1] - (data.loglam)[0]) * alog(10d) 
     sampl = 1
     SignalLin = uly_spect_alloc(DATA=flux, START=start, STEP=step, ERR=err, SAMP=1) 
  endif $
  else SignalLin = uly_spect_alloc(DATA=flux, WAVELEN=wave, ERR=err, GOODPIX=goodpix, SAMP=2)

  dof_factor = double(string(sxpar(h, 'DOF_FACT', COUNT=count)))
  if count eq 1 then SignalLin.dof_factor = dof_factor

  return, SignalLin
end

;==============================================================================
; reading routine of 1D and LSS spectra (format=3)
function uly_spect_read_lss, data, h, DISP_AXIS=disp_axis, QUIET=quiet

SignalLin = uly_spect_alloc(DATA=data, HEAD=h)
naxis = sxpar(h, 'NAXIS')


; ----------------------------------------------------------------------------
; search the dispersion axis, and the dispersion type (lin or log)
disp_type = -1
if n_elements(disp_axis) le 0 then disp_axis = 0
vacuum = 0  ; Set to 1 if the wavelengths are in vacuum

if disp_axis le 0 then begin ;  is there a standard WCS? (Vacuum wavelength)
    ctype = sxpar(h,'CTYPE*',COUNT=cnt)
    if cnt gt 0 then begin
        if ctype[0] eq 'WAVE-WAV' then begin ;if it is 2d we need the first CTYPE
            creator = sxpar(h,'CREATOR', COUNT=count) ; try to patch a Pleinpot bug
            if count gt 0 and creator eq 'Pleinpot     1' then ctype = 'AWAV'
;  in Elodie lib. 'WAVE' means air (in fits standarts 'AWAV' is air wavelength)
        endif
        disp_axis = 1 + where(strtrim(ctype,2) eq 'WAVE' or strmid(ctype,0,5) eq 'WAVE-') 
        if disp_axis gt 0 then begin
            if not keyword_set(quiet) then $
              message, /INFO, 'Dispersion axis is '+strtrim(disp_axis[0],2)+' (Vacuum wavelength)'
            vacuum = 1
        endif
    endif
endif

if disp_axis le 0 then begin ;  is there a standard WCS? (Air wavelength)
    disp_axis = 1 + where((strmid(ctype,0,4) eq 'AWAV') eq 1) 
    if disp_axis gt 0 then begin
        if not keyword_set(quiet) then $
          message, /INFO, 'Dispersion axis is '+strtrim(disp_axis[0],2)+' (Air wavelength)'
        if max(strmid(ctype[disp_axis-1],5,3) eq ['   ', 'WAV']) then begin
            disp_type = 0
            if not keyword_set(quiet) then message, /INFO, 'Dispersion axis is linear'
        endif else if strmid(ctype[disp_axis-1],5,3) eq 'LOG' then begin
            disp_type = 1
            if not keyword_set(quiet) then $
              message, /INFO, 'Dispersion axis is logarithmic'
        endif
    endif
endif

if disp_axis eq 0 then begin
    disp_axis = 1
    if naxis gt 1 and not keyword_set(quiet) then $
      message, 'Assume that the dispersion axis is 1 (X)'+$
               '(use the keyword DISP_AXIS to set it differently)', /INFO
endif else $
  if disp_axis eq 2 then *SignalLin.data = transpose(*SignalLin.data)

ax = strtrim(disp_axis[0],2)
crval = double(string(sxpar(h,'CRVAL'+ax)))

if disp_type lt 0 then begin    ; search the dispersion mode
    disp_type = 0
    if crval lt 10 then begin
        disp_type = 1
        if crval lt 5 then begin
            if not keyword_set(quiet) then $
              print, 'Assume that the dispersion is in log10',crval,'CRVAL'+ax
        endif else if not keyword_set(quiet) then $
          message, /INFO, 'Assume that the dispersion is in log (air wavelength)'
    endif else if not keyword_set(quiet) then $
      message, /INFO, 'Assume that the dispersion is linear (air wavelength)'
endif 

if disp_type ne 0 and disp_type ne 1 then $
  message, 'Do not handle this sampling (yet)'

; ----------------------------------------------------------------------------
; decode the spectral WCS
cdelt = double(string(sxpar(h,'CD'+ax+'_'+ax, COUNT=count)))
if count eq 0 then cdelt = double(string(sxpar(h,'CDELT'+ax)))

crpix = double(sxpar(h, 'CRPIX'+ax, COUNT=count))
if count eq 0 then crpix = 1d
crval = crval - (crpix - 1d) * cdelt ; wavelength of the 1st pixel

if (cdelt le 0.) or (crval le 0.) then $
  message,'WCS of the observations not set correctly'

if disp_type eq 1 and crval lt 5 then begin  
;   normally axis should be logn, but sometime it is log 10 ...
    if not keyword_set(quiet) then $
      print, 'Convert axis scale from log10 to log'
    crval *= alog(10d)
    cdelt *= alog(10d)
endif

if vacuum eq 1 then begin
    if not keyword_set(quiet) then $
      print, 'Wavelength in VACUUM ... approximately converted'
    if disp_type eq 0 then begin 
        crval /= 1.00028d 
        cdelt /= 1.00028d 
    endif else if disp_type eq 1 then crval -= 0.00028D
endif

; ----------------------------------------------------------------------------
; load the output structure

SignalLin.sampling = disp_type ; sampling in wavelength: lin/log
SignalLin.start = crval
SignalLin.step = cdelt

dof_factor = double(string(sxpar(h, 'DOF_FACT', COUNT=count)))
if count eq 1 then SignalLin.dof_factor = dof_factor

*SignalLin.hdr = h[3:*]      ; initialize the header
sxdelpar, *SignalLin.hdr, $     ; remove wcs and array specific keywords
  ['NAXIS1', 'CRVAL1', 'CRPIX1', 'CD1_1', 'CDELT1', 'CTYPE1', 'CROTA1', $
   'CD2_1', 'CD1_2', 'DATAMIN', 'DATAMAX', 'CHECKSUM' $
  ]

return, SignalLin

end

;============================================================================
function uly_spect_read, file_in, lmin, lmax,                           $
                         VELSCALE=velscale, SG=sg,                      $
                         ERR_SP=err_sp, SNR_SP=snr_sp, MSK_SP=msk_sp,   $
                         DISP_AXIS=disp_axis, FORMAT=format, QUIET=quiet

;; read the first extension that we hope contains a spectrum
;fits_read, file_in, data, h

; test if the file or unit argument is valid
if size(file_in, /TYPE) ne 3 and size(file_in, /TYPE) ne 7 then begin
    print, 'usage: ULY_SPECT_READ <filename>, ...'
    print, 'first argument must be a file name or unit number'
    return, 0
endif
file_inl = file_in  ; local copy of the argument
if size(file_in, /TYPE) eq 7 then begin
    if file_test(file_in) ne 1 then begin
        file_inl += '.fits'
        if file_test(file_inl) ne 1 then begin
            print, 'usage: ULY_SPECT_READ <filename>, ...'
            print, 'Error, file does not exist (' + file_in + ')'
            return, 0
        endif
    endif
endif

if n_elements(sg) eq 1 then begin
    if abs(sg) gt 10 then message, /INFO, $
      'Note that the SG (redshift) has an odd value: '+strtrim(sg,2)+$
      ' is it correct? (it should be a "z", not a "cz")'
endif

; read the first non-empty extension that we hope contains a spectrum
naxis = 0
nhdu = 0
status = 0
while naxis eq 0 and status eq 0 do begin       ; skip leading empty HDUs
    data = mrdfits(file_inl, nhdu, h, /SILENT, STATUS=status)
    if n_elements(h) eq 0 then begin
        print, 'Cannot access to the data in the file (invalid format?)'
        return, 0
    endif
    naxis = sxpar(h,'NAXIS')
    nhdu++
endwhile
if status ne 0 then begin
    print, 'Could not find a valid HDU in ', file_inl
    return, 0
endif

; switch to the appropriate reading routine
if n_elements(format) eq 0 then fmt = uly_spect_filefmt(h) $
else fmt = format
case fmt of
    1 : spect = uly_spect_read_sdss(data, h, ERR_SP=err_sp, SNR_SP=snr_sp, QUIET=quiet)
    2 : spect = uly_spect_read_tbl(data, h, ERR_SP=err_sp, SNR_SP=snr_sp, QUIET=quiet)
    3 : spect = uly_spect_read_lss(data, h, DISP_AXIS=disp_axis, QUIET=quiet)
;    4 : spect = uly_spect_read_sdss_dr9(data, h, FILE=file_inl, ERR_SP=err_sp, SNR_SP=snr_sp, QUIET=quiet) 
    else : begin
        print,'Could not recognize the format of FITS file'
        return, uly_spect_alloc(TITLE=file_in)
    end
endcase

ntot = n_elements(*spect.data)

if ntot eq 0 then begin
    print, 'No data read from FITS file'
    return, spect
endif

spect.title = file_in

; Handle the case when error  (or signal to noise) is read from another file
;    assume WCS are the same for error & data spectra
if n_elements(err_sp) gt 0 then begin
    testfile = FILE_INFO(err_sp)    
    if testfile.exists eq 1 then begin
        fits_read, err_sp, *spect.err, h_err 
        if disp_axis eq 2 then err = transpose(err)
    endif else $
       message, 'File:' + err_sp + ' does not exsists...'
    for i=0,n_elements((*spect.err)[0,*])-1 do begin ;for 2D case
       nans = where(finite((*spect.err)[*,i]) eq 0, cnt, COMP=fin)
       if cnt ne 0 then $
          if n_elements(fin) ne 1 then $ ;check if there were finite values
             (*spect.err)[nans,i] = max((*spect.err)[fin,i]) else $
                (*spect.err)[nans,i] = max((*spect.err)[*,i-1]) ; otherwise take the previous maximum error
    endfor
endif else if keyword_set(snr_sp) then begin
    testfile = FILE_INFO(snr_sp)
    if testfile.exists eq 1 then begin
        fits_read, snr_sp, err
        if disp_axis eq 2 then err = transpose(err)
        *spect.err = *spect.data / err
    endif else $
      message, 'SNR spectrum file not valid'
    for i=0,n_elements((err)[0,*])-1 do begin ;for 2D case
       neg = where(err[*,i] le 0, c, COMPLEM=pos)
       if c gt 0 then begin
          err[neg] = 1
          em = 10 * max((*spect.err)[pos,i])
          (*spect.err)[neg,i] = em
          message, 'The SNR spectrum '+strtrim(string(i),2)+ ' has '+strtrim(string(c),2)+$
                   ' negative or null values. Their error is set to '+$
                   strtrim(string(em),2), /INFO
       endif 
       nans = where(finite((*spect.err)[*,i]) eq 0, cnt, COMP=fin)
       if cnt ne 0 then (*spect.err)[nans,i] = max((*spect.err)[fin,i])
    endfor    
 endif

if n_elements(msk_sp) gt 0 then begin
    message, 'Read MASK spectrum ... NOT YET IMPLEMENTED'
endif

; Apply the shift to restframe if required
if n_elements(sg) gt 0 then begin  ; shift to rest-frame
    z1 = 1d + sg
    case spect.sampling of
        0 : begin
            spect.start /= z1
            spect.step /= z1
        end
        1 : spect.start -= alog(z1)
        2 : *spect.wavelen /= z1
    endcase
endif 

; Determine wavelength range in signal spectrum and extract the required region
if n_elements(lmin) gt 0 or n_elements(lmax) gt 0 then begin
    wr = uly_spect_get(spect, /WAVERANGE)
    if n_elements(lmin) gt 0 then wr[0] = min(lmin)
    if n_elements(lmax) gt 0 then wr[1] = max(lmax)
    if n_elements(velscale) eq 1 then begin
        wr[0] *= 1D - velscale/299792.458D/2D
        wr[1] *= 1D + velscale/299792.458D/2D
    endif
    spect = uly_spect_extract(spect, WAVERANGE=wr, /OVERWRITE)
    npix = (size(*spect.data))[1]
    spect.start = wr[0]
    spect.step = (wr[1]-wr[0])/float(npix)
endif

ntot = n_elements(*spect.data)  

; Check and eventually replace (with 0) the NaNs, and put in goodpixels
;   (we must be sure there is no NaNs before rebinning)
good = where(finite(*spect.data), cnt, COMPLEM=nans, NCOMPLEM=nnans)
if nnans gt 0 then begin
    if not keyword_set(quiet) then $
      print, 'The input spectrum contains'+string(n_elements(nans))+' NaNs ...'
    if n_elements(*spect.goodpix) eq 0 then *spect.goodpix = good $
    else begin
       maskI = bytarr(ntot)
       maskI[*spect.goodpix] = 1
       maskI[nans] = 0
       *spect.goodpix = where(maskI eq 1)
    endelse

;   patch 1 pix of the nan regions by replicating the edge value
;   in order to reduce the oscillations of the spline interpolation
    next = nans+1
    if next[n_elements(next)-1] eq n_elements(*spect.data) then $
      next[n_elements(next)-1] = nans[n_elements(next)-1]
 ;we connect here the last and the firs pixels in 2D data (not good)
    (*spect.data)[nans] = (*spect.data)[next]
    nans = where(finite(*spect.data) eq 0, cnt)
    if cnt gt 0 then begin
        prev = nans - 1
        if prev[0] lt 0 then prev[0] = nans[1]
        (*spect.data)[nans] = (*spect.data)[prev]
    endif
    nans = where(finite(*spect.data) eq 0, cnt)

    if cnt gt 0 then (*spect.data)[nans] = 0
endif

; do the same with the error, if exists
if n_elements(*spect.err) gt 0 then begin
    good = where(finite(*spect.err), cnt, COMPLEM=nans)
    if cnt eq 0 then begin
        if not keyword_set(quiet) then begin
            message, /INFO, 'The error spectrum does not contain finite values'+$
                     '... ignore it (ie. do as if no errors were given)'
        endif
        undefine, *spect.err
    endif else if cnt lt n_elements(*spect.err) then begin
        if not keyword_set(quiet) then $
          message, /INFO, 'The input spectrum contains'+string(n_elements(nans))+' NaNs ...'
        if n_elements(*spect.goodpix) eq 0 then *spect.goodpix = good $
        else begin
           maskI = bytarr(ntot)
           maskI[*spect.goodpix] = 1
           maskI[nans] = 0
           *spect.goodpix = where(maskI eq 1, cnt)
           if cnt eq 0 then begin
               if not keyword_set(quiet) then message, /INFO, $
                'No good pixels would be left after applying the MASK ... ignore the mask'
               undefine, *spect.goodpix
           endif
        endelse
        
;   patch 1 pix of the nan regions by replicating the edge value
;   in order to reduce the oscillations of the spline interpolation
        next = nans + 1
        if next[n_elements(next)-1] eq n_elements(*spect.err) then begin
            if n_elements(nans) gt 1 then begin
                nans = nans[0:n_elements(nans)-2]
                next = next[0:n_elements(nans)-1]
            endif else next = nans
        endif 
        (*spect.err)[nans] = (*spect.err)[next]

        nans = where(finite(*spect.err) eq 0, cnt)
        if cnt gt 0 then begin
            prev = nans - 1
            if prev[0] lt 0 then begin
               nans = nans[1:*]
               prev = prev[1:*]
            endif
            (*spect.err)[nans] = (*spect.err)[prev]
        endif
        nans = where(finite(*spect.err) eq 0, cnt)
        
        if cnt gt 0 then (*spect.err)[nans] = 0
    endif
endif

; rebin in log of wavelength
if n_elements(velscale) ne 0 then begin
    if n_elements(lmin) eq 0 then $
      spect = uly_spect_logrebin(spect, velscale, /OVER) $
    else $
      spect = uly_spect_logrebin(spect, velscale, WAVERANGE=lmin[0], /OVER)
endif

;after rebining different number of pxs
ntot = n_elements(*spect.data)
dim = size(*spect.data, /DIM)

; select goodpixels the good pixels are between lmin[i] and lmax[i] 
npix = (size(*spect.data))[1]
Pix_gal = spect.start + lindgen(npix) * spect.step
;Pix_gal = range(spect.start, spect.start+(npix-1d)*spect.step,npix)

if n_elements(lmin) eq 0 then lmn = [spect.start] else begin
    lmn = lmin
    if spect.sampling eq 1 then lmn = alog(lmn) 
endelse
if n_elements(lmax) eq 0 then lmx = [spect.start+spect.step*(npix-1)] else begin
    lmx = lmax
    if spect.sampling eq 1 then lmx = alog(lmx) 
endelse

good = 0L
for i = 0, n_elements(lmn) - 1 do begin
    good = [good, $
            where((Pix_gal gt lmn[i]) and (Pix_gal lt lmx[i]))]
endfor 

if (n_elements(*spect.goodpix) eq 0) and (n_elements(good) gt 1) then begin
   maskI = bytarr(ntot)
   maskI = reform(maskI, dim)
   maskI[good[1:*], *] = 1   
   *spect.goodpix = where(maskI eq 1)
endif else begin                ; have to combine with the previous mask
   maskI = bytarr(ntot)
   maskI[*spect.goodpix] = 1
   maskI = reform(maskI, dim)
   maskI[good[1:*],*] += 1
   *spect.goodpix = where(maskI eq 2)
endelse

sxaddpar, *spect.hdr, 'HISTORY', 'uly_spect_read, '''+strtrim(file_in,2)+'''

return, spect

end

;== end =======================================================================
