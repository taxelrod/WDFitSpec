;----------------------------------------------------------------------------
;+
; NAME:        
;              ULY_LINE
;
; PURPOSE:     
;              Define a spectral-line component 
;
; USAGE:
;              cmp = uly_line(wave, sigma, LSF=lsf_file, NAME=name, 
;              LL=lim_wave, SL=lim_sigma, WL=lim_weight)
;
; ARGUMENT:
;   WAVE:      Input. [Angstroem] Wavelength of the line
;   SIGMA:     Input. [km/s] Broadening of the line
;
; KEYWORDS:
;   LSF:       Name of the file containing a relative LSF to
;              be injected in the template. 
;
;   NAME:      A string, name of this component. By default a 
;              unique name is generated.
;
;   LL:        Permitted limits for the fitted wavelength (in Angstroem). 
;              If given, it must be an array of two values.
;
;   SL:        Permitted limits for the velocity dispersion, sigma (in km/s).
;              If given, it must be an array of two values.
;
;   WL:        Permitted limits for the weight of each model component
;              [min_w, max_w]; given not normalized.
;              If a single value is passed, take it for 'min_w'.
;              Note that this constraints works only in the case of multiple
;              components.
;
; DESCRIPTION:
;     Return a component to fit a single line (emission or absorption).
;     The model is either an unresolved line (Dirac, if SIGMA is not given)
;     or a Gaussian of standard deviation SIGMA, in km/s. (note that
;     the line may be convolved with a LOSVD in the fitting process
;     if ULYSS is called with KMOMENT > 1).
;     
;     The free parameters are the wavelength of the line, its gaussian
;     broadening and the weight of the component, determining the
;     height or depth (intensity) of the line.
;     In the future we may modify this procedure to make the 
;     Gauss-Hermite coefficients (h3 h4) of the line other free parameters.
;
;     This cmp can be passed to ULYSS.
;
;     The functions ULY_LINE_INIT and ULY_LINE_EVAL attached in this
;     file are automatically called by the fitting routines.
;
; OUTPUT:     
;     cmp structure or array of structures, for more detailed explanation
;     check ULY_FIT.
;
; EXAMPLE:
;     The following example shows how to fit the H & K lines in a
;     spectrum. Note that using an additive term (AD=0) is required
;     to fit the continuum.
;              spect = uly_spect_read(uly_root+'/data/m67.fits')
;              cmp1 = uly_line(3969e, 400e)
;              cmp2 = uly_line(3934e, 400e)
;              cmp = [cmp1, cmp2]
;              ulyss, spect, cmp, AD=0, /PLO
;
; SEE ALSO:
;     Definition of other components: ULY_SSP, ULY_TGM, ULY_STAR
;
; REQUIRED FUNCTION:
;     ULY_CMP_NAME ULY_SPECT_ALLOC ULY_SPECT_FREE ULY_SPECT_LSFCONVOL
;     ULY_SPECT_LOSVDCONVOL
;
; HISTORY:
;              Mina Koleva, 2008/05 created
;
;-
; CATEGORY:    ULY_UTIL
;------------------------------------------------------------------------------
function uly_line_eval, eval_data, pars

compile_opt idl2
on_error, 2

c = 299792.458d

pos = (alog(pars[0]) - eval_data.start) / eval_data.step

if pos lt 0 or pos ge eval_data.npix then return, fltarr(eval_data.npix)

;;;;;
;xborder = dindgen(eval_data.npix+1) - 0.5
;sigma_pix = pars[1] / (eval_data.step * c)
;y = (xborder -pos) / sigma_pix
;cumul = 0.5*(1+erf(y/sqrt(2d)))
;data0 = (shift(cumul,-1)-cumul[0:eval_data.npix-1])  
;;;;;

sigma_pix = pars[1] / (eval_data.step * c)
nline = ceil(10 * sigma_pix) 
xborder = round(pos) + dindgen(2*nline+2) - nline - 0.5d
cumul = 0.5 * (1 + erf((xborder-pos)/sigma_pix /sqrt(2d)))
data = fltarr(eval_data.npix)
p1 = round(pos-nline)
p2 = 0
if p1 lt 0 then begin
    p1 = 0
    p2 = -p1
endif
n = min([n_elements(cumul)-1,eval_data.npix-p1+p2])-1
data[p1] = (shift(cumul,-1)-cumul[0:2*nline])[p2:n]


; LSF injection (this is not the right way to do! TO BE CHANGED)
; (we should simply read the lsf file to find the value of sigma
;  at the position of the line and combine this value with sigma_pix above)
if n_elements(*eval_data.lsf_file) gt 0 then begin
   spect = uly_spect_alloc(DATA=data, START=eval_data.start, STEP=eval_data.step, SAMPLING=1)
   uly_spect_lsfconvol, *eval_data.lsf_file, spect
   data = *spect.data
   uly_spect_free, spect
endif

return, data

end

;------------------------------------------------------------------------------
function uly_line_init, cmp, WAVERANGE=lamrange, VELSCALE=velscale, QUIET=quiet

compile_opt idl2
on_error, 2

c = 299792.498d

cmp.eval_fun = 'uly_line_eval'

cmp.start = alog(double(lamrange[0]))
cmp.step = double(velscale)/c
cmp.npix = round((alog(double(lamrange[1])) - cmp.start) / cmp.step) + 1
cmp.sampling = 1s
(*cmp.para)[0].step = *(*cmp.para)[0].guess * cmp.step * 1.1

init_data = *cmp.init_data

eval_data = {lsf_file:ptr_new(*init_data.lsf_file),         $
            start:cmp.start, step:cmp.step, npix:cmp.npix}

cmp.eval_data = ptr_new(eval_data)

if n_elements(*(*cmp.para)[1].guess) ne 0 then begin ; guess sigma specified
   sigma_pix = *(*cmp.para)[1].guess / (cmp.step * c)
   if (*cmp.para)[1].limited[0] eq 0 then $
      if (sigma_pix gt 0.3) ne 0 then begin
      (*cmp.para)[1].limited[0] = 1
      (*cmp.para)[1].limits[0] = 0.3 * cmp.step * c; 0.3 pixels
   endif else $
      if *(*cmp.para)[1].guess lt  (*cmp.para)[1].limits[0] then $
         message, 'ULY_LINE_INIT: Failed for '+cmp.name+ $
                  ' because sigma is not within the limits'
   if (*cmp.para)[1].limited[1] eq 1 and $
      *(*cmp.para)[1].guess gt  (*cmp.para)[1].limits[1] then $
         message, 'ULY_LINE_INIT: Failed for '+cmp.name+ $
                  ' because sigma is not within the limits'      
   if sigma_pix le 0.3 then (*cmp.para)[1].fixed = 1S
endif else begin
   if (*cmp.para)[1].limited[0] eq 0 then begin
      (*cmp.para)[1].limited[0] = 1
      (*cmp.para)[1].limits[0] = 0.3 * cmp.step * c; 0.3 pixels
      if (*cmp.para)[1].limited[1] eq 0 then  $
         *(*cmp.para)[1].guess = cmp.step * c  $
      else $
         *(*cmp.para)[1].guess = min([cmp.step*c, (*cmp.para)[1].limits[1]])
   endif
endelse

if *(*cmp.para)[1].guess / (cmp.step * c) ge 0.3 then (*cmp.para)[0].step /= 4.

(*cmp.para)[1].step = velscale/5.

return, cmp  

end

;------------------------------------------------------------------------------
function uly_line, wave, sigma,  $
                   LSF=lsf_file, $
                   NAME=name,    $
                   LL=lim_wave,  $
                   SL=lim_sigm,  $
                   WL=lim_weight

compile_opt idl2
on_error, 2

; Note
; This component is fitting a single emission line represented by a
; Dirac function, eventually convolved with a Gaussian of standard
; deviation SIGMA.
; If SIGMA=0 the precision of the line position is limited to 1 pix:
; this is an important limitation of the program ... ULYSS is not good
; for fitting unresolved lines.
; The line may be also shifted and broadened in ULY_FIT_LIN by the
; LOSVD convolution in the same time as all the other components.

if n_elements(wave) eq 0 then begin
    print, 'Usage: ULY_LINE, wave, ...'
    message, 'Error, at minimum one wavelength must be passed', /INFO
    return, 0
endif

s_wl = size(lim_weight)

namec = uly_cmp_name(name)

init_data = {lsf_file:ptr_new(lsf_file)}

if s_wl[0] gt 0 then if not (s_wl[0] eq 1 and s_wl[1] eq 2) then $ 
  message, 'The weight limits have to be of the type arr(2)'
    
case n_elements(lim_weight) of
   1: lw = [double(lim_weight), (machar(/DOUBLE)).xmax]
   2: lw = double(lim_weight)
   else : lw = [-(machar(/DOUBLE)).xmax, (machar(/DOUBLE)).xmax]
endcase

cmp = {name:namec, $
       descr:'emission line', $
       init_fun:'uly_line_init', $
       init_data:ptr_new(init_data), $
       eval_fun:'uly_line_eval', $
       eval_data:ptr_new(), $
       para:ptr_new(/ALLO), $
       start:0d, $
       step:0d, $
       npix: 0l, $
       sampling:-1s, $
       mask:ptr_new(), $
       weight:0d, $
       e_weight:0d, $
       l_weight:0d, $
       lim_weig:lw $
      }

; A single free parameter: The wavelength of the line
dwave = double(wave)
*cmp.para = [{name:'lambda', unit:'Angstrom', guess:ptr_new(double(dwave)), $
              step:1D, limits:[0d,0d], limited:[0S,0S], fixed:0S,           $
              value:0D, error:0D, dispf:''},                                $
             {name:'sigma', unit:'km/s', guess:ptr_new(double(sigma)),      $
              step:1D, limits:[0d,0d], limited:[0S,0S], fixed:0S,           $
              value:0D, error:0D, dispf:''}]

if n_elements(lim_wave) gt 0 then begin
    (*cmp.para)[0].limited = [1S, 1S]
    (*cmp.para)[0].limits = lim_wave
endif
if n_elements(lim_sigm) gt 0 then begin
    (*cmp.para)[1].limited = [1S, 1S]
    (*cmp.para)[1].limits = lim_sigm
endif

; return the structure describing this component
return, cmp

end


