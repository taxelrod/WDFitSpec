;+
; NAME:
;                 ULY_TGM
; PURPOSE:
;                 Define a TGM component
; USAGE:
;                 cmp = uly_tgm(MODEL_FILE=model_file, LSF=lsf_file,      
;                               TG=t_guess, LG=l_guess, ZG=z_guess,   
;                               TL=t_limits, LL=l_limits, ZL=z_limits,
;                               NAME=name, FIXPAR=fixpar, WL=lim_weight)
;
; KEYWORDS:
;   MODEL_FILE    Name of the FITS file containing the TGM model, by default 
;                 use: /models/elodie32_flux_tgm.fits.
;
;   LSF           Name of the file containing a relative LSF to be injected in 
;                 the template. 
;
;   TG            [K] Guess value for the effective temperature (Teff), by
;                 default = 7000.0. If an array is passed, the analysis will 
;                 be repeated starting from each element to find the optimal 
;                 solution.
;   
;   LG            [cm/s^2] Guess value for the surface gravity (Logg), by 
;                 default = 2.0, it can also be given as an array to search 
;                 the absolute minimum.
;
;   ZG            [dex] Guess value for the metallicity ([Fe/H]), by 
;                 default = 0.0, it can also be given as an array
;                 to search the best minimum from the grid of guesses.
;
;   TL            [K], dblarr(2), limits for the Teff,
;                 by default = [3100.0,40000.0] see  ULY_TGM_INIT.
;
;   LL            [dex],dblarr(2), limits for the Logg,
;                 by default = [-0.25, 5.9] see ULY_TGM_INIT.
;
;   ZL            [dex],dblarr(2), limits for the [Fe/H],
;                 by default = [-2.5, 1.0] see ULY_TGM_INIT.
;
;   NAME          A string, name of this component. By default a unique name is 
;                 generated.
;   
;   FIXPAR        An array used to fix some the atmosphere parameters,
;                 0 means that the parameter is free, and 1 that it is fixed.
;                 A fixed parameter keeps its initial value.
;                 The parameters are specified in the following order:
;                 (Teff, Logg, [Fe/H]).
;          
;   WL            Limits for the weight of the component [min_w, max_w].
;                 The weight is in data units over component units.
;                 By default the weight is constrained to be positive.
;                 This constraint is applicable when multiple components are
;                 fitted (composite spectrum), it is ignored for a single
;                 component.
;
;   /REBIN_COEF   Computation mode, see the description below.
;       
; DESCRIPTION:
;     Return a TGM structure describing a TGM component that can be 
;     passed to program ULYSS. 
;
;     TGM stands for Temperature-Gravity-Metallicity: A TGM component 
;     specifies an interpolated stellar spectrum evaluated as a function 
;     of its atmospheric parameters (Teff, Logg and [Fe/H]) that was 
;     interpolated over a stellar spectral library (e.g. ELODIE 3.2) 
;     or over a grid of synthetic spectra. 
;
;     The file passed to ULY_TGM contains the coefficients of polynomials 
;     depending on the atmospheric parameter. Each wavelengh bin is represented
;     separately, therefore the coefficients for each term of the development
;     form a spectrum. Depending on the interpolator, there are 14 to about
;     25 terms. The first polynomial interpolator was presented in 
;     Prugniel and Soubiran 2001, A&A 369, 1048 for the ELODIE spectral
;     libraries, it was updated in subsequent papers, in particular in 
;     Wu et al. 2011, A&A 525, 71. Interpolators for the MILES library were
;     presented in Prugniel et al. 2011, A&A 531, 169, and Sharma et
;     al. 2015, submitted. Other interpolators are available on the
;     ULySS site.
;
;     The interpolated spectra are evaluated using the interpolator, and
;     then rebinned to the WCS of the spectrum to analyse. 
;
;     REBIN_COEF option:
;     Rather than rebinning the interpolated spectrum, it first rebin the 
;     coefficients, and then compute the interpolated spectrum. This option is 
;     slightly less precise (because of the spline interpolation of the several
;     coefficients), but is in some cases faster. This is the case when (i)
;     a series of several spectra are fitted with the same CMP, and (ii) when
;     the spectrum to analyse has a sampling much coarser than the interpolator.
;     (Note that the impact on the precision of the interpolated spectrum
;     is generally fully acceptable; it was the only available option before
;     the release 1.3 of ULySS.)
;     
; OUTPUT:
;     TGM cmp struct, for detailed explanation please check uly_fit.pro        
;
; REQUIRED FUNCTION:
;                 ULY_CMP_NAME
;                         
; EXAMPLE:
;     Fit a CFLIB star: 
;     cmp = uly_tgm()
;     star = uly_root+'/data/cflib_114642.fits'
;     ulyss, star, cmp, lmin=3900.0, lmax=6800.0, /PLOT
;
; HISTORY:
;                 Creation Yue WU 2008/06/18
;
;-
; CATEGORY:       ULY_TGM
;------------------------------------------------------------------------------
function uly_tgm, MODEL_FILE=model_file,                  $
                  LSF=lsf_file,                           $
                  TL=t_limits, LL=l_limits, ZL=z_limits,  $
                  TG=t_guess, LG=l_guess, ZG=z_guess,     $
                  NAME=name, FIXPAR=fixpar, WL=lim_weight,$
                  REBIN_COEF=rebin_coef
                                                                    
compile_opt idl2
on_error, 2

common uly_path, uly_root

if n_elements(model_file) eq 0 then $
  model_file = uly_root + '/models/elodie32_flux_tgm.fits'

if size(model_file, /TYPE) ne 7 then begin
    message, 'Argument model_file must be a filename', /INFO
    return, 0
endif
file = strtrim(model_file, 2)
if file_test(file) ne 1 then begin
    file += '.fits'
    if file_test(file) ne 1 then begin
        message, 'Argument model_file must be a filename ('+model_file+')', /INFO
        return, 0
    endif
endif

if n_elements(rebin_coef) eq 0 then rebin_coef = 0

if n_elements(t_guess) eq 0  then t_guess = 7000.0d
if n_elements(l_guess) eq 0  then l_guess = 2.0d
if n_elements(z_guess) eq 0  then z_guess = 0.0d

s_t = size(t_limits)
if s_t[0] gt 0 then $
  if not s_t[0] eq 2 then $
  message, 'The Teff limits have to be of the type arr(2)'

s_l = size(l_limits)
if s_l[0] gt 0 then $
  if not s_l[0] eq 2 then $
  message, 'The Logg limits have to be of the type arr(2)'

s_z = size(z_limits)
if s_z[0] gt 0 then $
  if not s_z[0] eq 2 then $
  message, 'The metallicity limits have to be of the type arr(2)'

s_wl = size(lim_weight)
if s_wl[0] gt 0 then $
   if not s_wl[0] eq 2 then $
   message, 'The weight limits have to be of the type arr(2)'

; create the struct describing TGM component
init_data = {model:file, lsf_file:ptr_new(lsf_file), rebin_coef:rebin_coef}
descr = ''
if n_elements(file) eq 1 then descr += 'model:' + file + ' ' 
if n_elements(lsf_file) eq 1 then descr += 'lsf:' + lsf_file + ' '
if rebin_coef eq 1 then descr += 'rebin_coef ' 

namec = uly_cmp_name(name)

cmp = {name:namec, $
       descr:descr, $
       init_fun:'uly_tgm_init', $
       init_data:ptr_new(init_data), $
       eval_fun:'', $
       eval_data:ptr_new(), $
       para:ptr_new(/ALLO), $
       start:0d, $
       step: 0d, $
       npix: 0l, $
       sampling:-1s, $
       mask:ptr_new(), $
       weight:0d, $
       e_weight:0d, $
       l_weight:0d, $
       lim_weig:(machar(/DOUBLE)).xmax*[0,1] $
      }

; load the parameters and data in the component struct
*cmp.para = replicate({name:'', unit:'', guess:ptr_new(), step:1D-2, $
                       limits:[0d,0d], limited:[1,1], fixed:0S, $
                       value:0D, error:0D, dispf:''}, $
                      3)       

(*cmp.para)[0].name = 'Teff'
(*cmp.para)[0].unit = 'K'
(*cmp.para)[0].dispf = 'exp'  ; says to display exp(para)
(*cmp.para)[1].name = 'Logg'
(*cmp.para)[1].unit = 'cm/s2'
(*cmp.para)[2].name = 'Fe/H'
(*cmp.para)[2].unit = 'dex'
 
if n_elements(t_limits) gt 0 then (*cmp.para)[0].limits = alog(t_limits)
if n_elements(l_limits) gt 0 then (*cmp.para)[1].limits = l_limits
if n_elements(z_limits) gt 0 then (*cmp.para)[2].limits = z_limits

for k=0,2 do (*cmp.para)[k].guess = ptr_new((*cmp.para)[k].limits[0])

if n_elements(t_guess) ge 1 then $
  *(*cmp.para)[0].guess = alog(double(string(t_guess)))
if n_elements(l_guess) ge 1 then $
  *(*cmp.para)[1].guess = double(string(l_guess))
if n_elements(z_guess) ge 1 then $
  *(*cmp.para)[2].guess = double(string(z_guess))

(*cmp.para).step = 0.01d  
(*cmp.para)[0].step = 0.005d  ; derivation step in Teff = 0.5 %

if n_elements(fixpar) gt 0 then $
  (*cmp.para)[0:n_elements(fixpar)-1].fixed = fixpar

if n_elements(lim_weight) gt 0 then $
  cmp.lim_weig = double(lim_weight)

return, cmp

end
