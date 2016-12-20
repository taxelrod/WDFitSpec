;+
; NAME:            ULY_TGM_EVAL
; 
; PURPOSE:           
;                  Evaluate a TGM model spectrum.
; 
; USAGE:
;                  array = uly_tgm_eval(eval_data, para)
; 
; ARGUMENTS:
;   EVAL_DATA:     The interpolated stellar spectrum model coefficients, 
;                  it contains three dimensions of data, specifying to
;                  warm, hot and cold stellar regions respectively, in
;                  fact, they are the coefficients of the evaluating 
;                  polynomial function.
;
;   PARA:          Input 3 atmospheric parameters (log(Teff), Logg, [Fe/H]),
;                  they will be passed to ULY_TGM_MODEL_PARAM to
;                  compute out the parameters which will be applied
;                  to the polynomial function for evaluating
;                  the model spectrum at this given stellar
;                  atmospheric parameters (Teff, log g and [Fe/H]).
;
; DESCRIPTION:
;      Evaluate a TGM model spectrum which is interpolated as a
;      function of atmospheric parameters (effective temperature,
;      surface gravity and [Fe/H]). Each wavelength point is evaluated
;      by using polynomials.
;
;      This function is automatically called by the minimization routines:
;      ULY_FIT, ULY_FIT_LIN to evaluate a TGM model (see ULY_TGM).
; 
; OUTPUT:          
;                  Evaluated stellar model spectrum               
;
; REQUIRED FUNCTION: 
;                  ULY_TGM_MODEL_PARAM, compute the parameters 
;       
; AUTHOR:
;                  Yue WU, 2008/06/18
;
;-
; CATEGORY:        ULY_TGM
;------------------------------------------------------------------------------
function uly_tgm_model_param, version, teff, gravi, fehi  
                                                                 
compile_opt idl2, hidden
on_error, 2

if version eq 4 then param = dblarr(16) $
else message, 'Version'+string(version)+' of TGM not supported'

tt = teff*1d0
tt2 = tt*tt
tt3 = tt2*tt                            
grav = gravi * 1d0
grav2 = grav*grav
grav3 = grav2*grav
feh = fehi * 1d0

param[0]  = tt3*grav3                                                      
param[1]  = tt3*grav2
param[2]  = tt3*grav
param[3]  = tt3
param[4]  = tt2*grav3               
param[5]  = tt2*grav2                                                 
param[6]  = tt2*grav
param[7]  = tt2
param[8]  = tt*grav3                                                 
param[9]  = tt*grav2                                              
param[10] = tt*grav                                                 
param[11] = tt
param[12] = grav3
param[13] = grav2
param[14] = grav
param[15] = 1

return, param                                                    
end   

;------------------------------------------------------------------------------
function uly_tgm_eval, eval_data, para

spec_coef = eval_data.spec_coef
spec_npix = (size(spec_coef, /DIM))[0]

teff = para[0]
grav = para[1]

param = uly_tgm_model_param(eval_data.version, teff, grav, para[2])
np =  (size(param,/DIM))[0] - 1

tgm_model_evalhc = spec_coef[*,0:np] # param[*]

; Resample in wavelength if necessary (when option REBIN_COEF is not given)
if eval_data.sampling ne eval_data.mod_samp or $
   eval_data.start ne eval_data.mod_start or $
   eval_data.step ne eval_data.mod_step or $
   eval_data.npix ne spec_npix then begin
   spec = uly_spect_alloc(DATA=tgm_model_evalhc, START=eval_data.mod_start, $
                          STEP=eval_data.mod_step, SAMPLING=eval_data.mod_samp)
   wrange = eval_data.start + [0d, eval_data.npix * eval_data.step]
   if eval_data.sampling eq 1 then wrange = exp(wrange)
   c = 299792.458d              ; Speed of light in km/s
   velscale = eval_data.step * c
   if eval_data.sampling eq 1 then $
      spec = uly_spect_logrebin(spec, velscale, WAVERANGE=wrange, /OVER) $
   else $
      spec = uly_spect_linrebin(spec, eval_data.step, WAVERANGE=wrange, /OVER)
   tgm_model_evalhc = *spec.data
   uly_spect_free, spec
endif

if eval_data.calibration eq 'C' then begin ; multiply by a black-body spectrum
  n = n_elements(tgm_model_evalhc)
  if eval_data.sampling eq 1 then $
    wavelength = exp(eval_data.start + dindgen(n) * eval_data.step) $
  else $
    wavelength = eval_data.start + dindgen(n) * eval_data.step 
  w = 5550.
  c3 = 1.43883d8 / 5550d0 / exp(para[0])
  c1 = 3.74185d19 / 5550d0^5
  if c3 lt 50. then bbm = (c1 / (exp(c3)-1.)) else bbm = (c1 * exp(-c3))

  c3 = 1.43883d8 / wavelength / exp(para[0])
  c1 = 3.74185d19 / wavelength^5 / bbm

  n1 = where(c3 lt 50, cnt, COMPLEM=n2)
  if cnt gt 0 then tgm_model_evalhc[n1] *= (c1[n1] / (exp(c3[n1])-1.))
  if cnt lt n then tgm_model_evalhc[n2] *= (c1[n2] * exp(-c3[n2]))

endif

; Convolve LOSVD in case giving lsf_file
; (if REBIN_COEF is given, the LSF injection shall be made at initialization.)
if n_elements(eval_data.lsf) gt 0 and eval_data.lsf ne 'no_lsf' then begin
    spec = uly_spect_alloc(DATA=tgm_model_evalhc, START=eval_data.start, $
                           STEP=eval_data.step, SAMPLING=1)
    uly_spect_lsfconvol, eval_data.lsf, spec
    tgm_model_evalhc = *spec.data
    uly_spect_free, spec
endif

return, tgm_model_evalhc

end

