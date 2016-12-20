;+
; NAME:
;              ULY_FIT_INIT
; PURPOSE:
;              Call the initialization function for each component of the fit
;              
; USAGE:
;              status = uly_fit_init(cmp, WAVERANGE=lamrange, 
;                                         VELSCALE=velscale, QUIET=quiet)
; ARGUMENTS:
;     cmp : The array describing the model to fit, see ULY_FIT.
;
; KEYWORDS:
;     WAVERANGE : (Input 2 elements double precision array) 
;                 Wavelength range, in Angstrom.
;     VELSCALE  : (input) Size of the pixel in km/s.
;     QUIET     : Set this keyword to suppress some informational messages.
;
; RETURN VALUE:
;     status is an error status, 0 for a sucessful completion.
;
; DESCRIPTION:
;     The model to fit is described in the cmp array that can be created
;     by, e.g., ULY_SSP. Each element of this array is a structure describing 
;     one component of the fit. One of the member of this structure is
;     the name of the initialization function, and ULY_FIT executes these
;     functions for each of the components.
;
;     The initialization consists in operations that can be made once at
;     the beginning to save computations at the time of the 'evaluation'.
;     For example, the initialization may read a file, resample a model
;     grid as the observation ...
;
;     ULY_FIT_INIT must be called after reading the observation (because
;     it needs to know the sampling and wavelength range) but before 
;     executing the fitting procedure ULY_FIT.
;
; HISTORY:
;              Philippe Prugniel, 2008/05 created
;-
; CATEGORY:    ULY
;------------------------------------------------------------------------------
function uly_fit_init, cmp, WAVERANGE=lamrange, VELSCALE=velscale, $
                       QUIET=quiet, NOREDO=noredo
                        
compile_opt idl2
on_error, 2
; The 'NOREDO' keyword is used to prevent recursive re-entry

wr = double(lamrange)
for i=0, n_elements(cmp)-1 do begin
    if strtrim(cmp[i].init_fun,2) eq '' then begin
        cmp[i].start = alog(wr[0])
        cmp[i].step = velscale / 299792.458d
        cmp[i].npix = 1 + round(alog(double(wr[1])/double(wr[0])) / cmp[i].step)
        cmp[i].sampling = 1
        *cmp[i].eval_data = cmp[i]
    endif else $
      cmp[i] = call_function(cmp[i].init_fun, cmp[i], WAVERANGE=wr, $
                             VELSCALE=velscale, QUIET=quiet)

    if cmp[i].start gt alog(wr[0]) then wr[0] = exp(cmp[i].start)
    if cmp[i].start+cmp[i].step*(cmp[i].npix-1) lt alog(wr[1]) then $
       wr[1] = exp(cmp[i].start+cmp[i].step*(cmp[i].npix-1))

    if cmp[i].sampling ne 1 then begin
        message, 'Some components of the model are not log-wave sampled ... ' + $
          strjoin(cmp[i].name,','), /CONT
        return, 1
    endif

endfor

start_min = min(cmp.start, MAX=start_max)
cmp_end = cmp.start + (cmp.npix-1)*cmp.step
end_min = min(cmp_end, MAX=end_max)

; If all the components do not have the same wavelength range, redo
; the initialization with the common range. This happens if the first
; CMP is not the one with the narrowest range.

; This is not very elegant, because we have to redo possibly heavy
; initialization computations, but it is the only way we have now.
; Hopefully this situation is not frequent, and it can be avoided by 
; selecting an appropriate range in the observation we want to analyse.

if start_max-start_min gt 0.02*velscale/299792.458d or $
  end_max-end_min gt 1.02*velscale/299792.458d then begin

    wr = exp([start_max, end_min])
    message, 'All the model components do not cover the same wavelength range...', /INFO

    if not keyword_set(noredo) then begin
        message, 'Initialization is redone with range ('+strjoin(strtrim(wr,2),',')+') ...', /INFO        
        return , uly_fit_init(cmp, WAVERANGE=wr, VELSCALE=velscale, QUIET=quiet, /NOREDO)
    endif

endif

return, 0

end
;--- end ----------------------------------------------------------------------
