;+
; NAME:
;               ULY_TGM_INIT
; PURPOSE: 
;               Initialize a TGM component
;             
; USAGE:
;               cmp_new = uly_tgm_init (cmp, WAVERANGE=lamrange, $
;                                       VELSCALE=velscale, QUIET=quiet)
;
; ARGUMENTS:   
;   CMP:        Component defined by using ULY_TGM, please check it for
;               more details.
;
; KEYWORDS:
;   WAVERANGE:  Wavelength range used when the model will be
;               evaluated with ULY_TGM_EVAL.
;
;   VELSCALE:   Size of one pixel in km/s. Used when the model will be
;               evaluated and log-rebinned with ULY_TGM_EVAL.
;
;   QUIET:      verbosity control
;
; DESCRIPTION: 
;      This initialization function is automatically executed by 
;      ULY_FIT_INIT when a component has been defined by using
;      ULY_TGM. A TGM model component can generate a spectrum at given 
;      stellar atmospheric parameters (Teff,log g and [Fe/H]). 
;
;      The name of this function is recorded as a member of
;      cmp.init_fun in the TGM component structure (see ULY_FIT for a
;      complete description of the component structure).
;
;      ULY_TGM_INIT consists in loading the basic TGM model
;      coefficients from a FITS file whose name is saved in
;      *cmp.init_data (by default
;      /models/elodie32_flux_tgm.fits), and also loads in the cmp 
;      structure the WCS.
;
; NOTE
; ULY_INIT_TGM is also called by ULY_TGM_EXTR where we use the keywords
; SAMPLING, WAVERANGE, STEP, which are more general because they can make
; either lin- or log-sampled spectra.
; This more general approach may be used in ULY_FIT_INIT in the future, when
; we will implement fitting in other than log- (not however that it is not
; fully general, because uneven sampling would not be supported).
;
; HISTORY:
;               Creation Yue WU 2008/06/18
; 
;-
; CATEGORY:     ULY_TGM
;---------------------------------------------------------------------


;==================================================================
; Load the tgm model's coefficients array
; Arguments:
;        filename, model file name
;
; This version of the code is ONLY for version 4 of the TGM files!
;
function uly_tgm_coef_load, filename

compile_opt idl2, hidden
on_error, 2

if n_elements(filename) eq 0 then begin
    message, 'no input model file name...', /INFO
    print, 'Usage: uly_tgm_model_load, <filename>, ...'
    return, 0
endif

fits_read, filename, data, hdr, MESSAGE=mess0  

if mess0 ne ''  then begin
    message, 'Failed to read FITS file...'+mess0, /INFO
    return, 0
endif

crval1 = double(string(sxpar(hdr, 'CRVAL1')))
cdelt1 = double(string(sxpar(hdr, 'CDELT1')))
naxis1 = string(sxpar(hdr, 'NAXIS1'))
ctype1   = strtrim(sxpar(hdr, 'CTYPE1', /SILENT, COUNT=cnt), 2)
if cnt eq 0 then ctype1 = 'AWAV'

if ctype1 eq 'AWAV' then samp = 0 $
else if ctype1 eq 'AWAV-LOG' then samp = 1 $
else message, 'Unsupported CTYPE1 value ' + ctype1

spec_coef = [[[data]]]
mask = bytarr(naxis1) + 1
goodpix = where (mask eq 1, cnt)
if cnt eq 0 then undefine, goodpix
s = uly_spect_alloc(START=crval1, STEP=cdelt1, SAMPLING=samp, DATA=spec_coef, GOODPIX=goodpix)

return, s
end


;====================================================
function uly_tgm_init, cmp, WAVERANGE=lamrange, VELSCALE=velscale, $
                            SAMPLING=sampling, STEP=step, QUIET=quiet
compile_opt idl2
on_error, 2

init_data = *cmp.init_data

cmp.eval_fun = 'uly_tgm_eval'

if size(init_data.model, /TYPE) eq 7 then begin
   if file_test(init_data.model) ne 1 then begin
      print, 'Error, model file does not exist!'
   endif
endif

fits_read, init_data.model, data, hdr, /HEADER_ONLY

uly_type = strtrim(sxpar(hdr, 'ULY_TYPE', /SILENT, COUNT=cnt), 2)
if cnt eq 1 and uly_type ne 'TGM' then $
   message,'Invalid model file, expect ULY_TYPE=TGM, get '+uly_type

naxis1   = string(sxpar(hdr, 'NAXIS1'))
naxis2   = sxpar(hdr, 'NAXIS2', /SILENT, COUNT=cnt)
if cnt eq 0 then $
   message,'Invalid model file, only one axis '
ctype1   = strtrim(sxpar(hdr, 'CTYPE1', /SILENT, COUNT=cnt), 2)
if cnt eq 0 then ctype1 = 'AWAV'
crval1   = double(string(sxpar(hdr, 'CRVAL1')))
cdelt1   = double(string(sxpar(hdr, 'CDELT1')))

version  = sxpar(hdr, 'INTRP_V', /SILENT, COUNT=cnt)
if cnt eq 0 then version = 1  ; there was no version KW at the beginning

calibration  = strtrim(sxpar(hdr, 'INTRP_C', /SILENT, COUNT=cnt), 2)

if naxis2 ne 16 then $
      message,'Invalid interpolator format for version='+strtrim(version,2)

if ptr_valid(cmp.eval_data) then begin
;  this cmp was yet initialized, we shall check whether we have to renitialize
   if n_elements(velscale) ne 0 and n_elements(lamrange) ne 0 then begin
      cmp_velsc = cmp.step * 299792.458d0
      if velscale eq cmp_velsc then begin
         cmp_npix = (size((*cmp.eval_data).spec_coef,/DIM))[0]
         cmp_wrang = cmp.start + [0,cmp_npix-1]*cmp.step
         if cmp.sampling eq 1 then cmp_wrang = exp(cmp_wrang)
         d = intarr(naxis1)
         s0 = uly_spect_alloc(START=crval1, STEP=cdelt1, SAMPLING=0, DATA=d)
         wr0 = uly_spect_get(s0, /WAVERANGE)
         wr = [max([wr0[0],lamrange[0]]), min([wr0[1],lamrange[1]])]
         s0 = uly_spect_logrebin(s0, velscale, WAVERANGE=wr, /OVER, /FLUX) 
         wr = uly_spect_get(s0, /WAVERANGE)
         uly_spect_free, s0
         if wr[0] eq cmp_wrang[0] and wr[1] eq cmp_wrang[1] then begin
            return, cmp
         endif
      endif
   endif
endif


; Set the parameters' limits
if total((*cmp.para)[0].limits ne [0d,0d]) eq 0 then begin
   teffrange = double(string(sxpar(hdr, 'TEFFLIM*', COUNT=cnt)))
   if cnt ne 2 then teffrange = [1000.0, 100000.]
   (*cmp.para)[0].limits = alog(teffrange) ; Teff limits in log
endif
if total((*cmp.para)[1].limits ne [0d,0d]) eq 0 then begin
   loggrange = double(string(sxpar(hdr, 'LOGGLIM*', COUNT=cnt)))
   if cnt ne 2 then loggrange = [6.0, 9.0]
   (*cmp.para)[1].limits = loggrange
endif
if total((*cmp.para)[2].limits ne [0d,0d]) eq 0 then begin
   fehrange = double(string(sxpar(hdr, 'FEHLIM*', COUNT=cnt)))
   if cnt ne 2 then fehrange = [-2.5, 1.0]
   (*cmp.para)[2].limits = fehrange
endif

s = uly_tgm_coef_load(init_data.model) 

wr = crval1 + [0, (naxis1-1)*cdelt1]
if ctype1 eq 'AWAV-LOG' then wr = exp(wr)
if n_elements(lamrange) gt 0 then wr[0] = max([lamrange[0],wr[0]])
if n_elements(lamrange) gt 1 then wr[1] = min([lamrange[1],wr[1]])

if init_data.rebin_coef eq 1 then begin
   sampl = 1
   if n_elements(sampling) eq 1 then sampl = sampling 
   if sampl eq 1 then begin
      undefine, velsc
      if n_elements(velscale) eq 1 then velsc = velscale $
      else if n_elements(step) eq 1 then velsc = step * 299792.458d 
      s = uly_spect_logrebin(s, velsc, WAVERANGE=wr, /OVER) 
   endif else if sampl eq 0 then begin
      s = uly_spect_linrebin(s, step, WAVERANGE=wr, /OVER) 
   endif
endif else begin
   s = uly_spect_extract(s, WAVERANGE=wr, /OVER) 
endelse

; Description of the model coefs stored in cmp (wavelength range, sampling and step)
mod_samp = s.sampling
mod_start = s.start
mod_step = s.step
spec_coef = *s.data

; Determine the WCS of the cmp itself
cmp.sampling = 1 ; by default the cmp will be log-sampled
if n_elements(sampling) eq 1 then cmp.sampling = sampling

if cmp.sampling eq mod_samp then cmp.step = mod_step 
if n_elements(velscale) eq 1 then begin
   if cmp.sampling ne 1 then message, 'Inconsistency in the arguments'
   cmp.step = velscale/299792.458d0
endif
if n_elements(step) eq 1 then cmp.step = step

resam = 0
if n_elements(lamrange) gt 0 then if wr[0] ne lamrange[0] then resam = 1
if n_elements(lamrange) gt 1 then if wr[1] ne lamrange[1] then resam = 1

if mod_samp ne cmp.sampling or mod_step ne cmp.step or resam eq 1 then begin
   s = uly_spect_extract(s, ONED=0, /OVER)
   if cmp.sampling eq 0 then begin
      undefine, step
      if cmp.step ne 0 then step= cmp.step
      s = uly_spect_linrebin(s, step, WAVERANGE=lamrange, /OVER)
   endif else if cmp.sampling eq 1 then begin
      undefine, velscale
      if cmp.step ne 0 then velscale = cmp.step * 299792.458d0
      s = uly_spect_logrebin(s, velscale, WAVERANGE=lamrange, /OVER)
   endif else message, 'Cannot yet resample to sampling=2'
endif

cmp.start = s.start
cmp.step = s.step
cmp.npix = (size(*s.data, /DIM))[0]
cmp.sampling = s.sampling

if n_elements(*s.goodpix) gt 0 then begin
   ptr_free, cmp.mask
   cmp.mask = ptr_new(bytarr(cmp.npix)) 
   (*cmp.mask)[*s.goodpix] = 1 
endif

uly_spect_free, s
ptr_free, cmp.eval_data

if n_elements(*init_data.lsf_file) gt 0 then lsf = *init_data.lsf_file $
else lsf = 'no_lsf'

cmp.eval_data = ptr_new({spec_coef:spec_coef, $
                         start:cmp.start, step:cmp.step, npix:cmp.npix, sampling:cmp.sampling, $
                         mod_samp:mod_samp, mod_start:mod_start, mod_step:mod_step, $
                         lsf:lsf, version:version, calibration:calibration})

return, cmp

end

