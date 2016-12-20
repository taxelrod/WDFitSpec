;+
; NAME:        ULY_STAR_INIT
;
; PURPOSE:     Initialize a single star fit component 
;
; USAGE:
;              cmp = uly_star_init(cmp, WAVERANGE=waverange, 
;                                  VELSCALE=velscale, /QUIET)
;
; ARGUMENTS:   
;   CMP:       A component defined using ULY_STAR, please check it for more
;              details.
;
; KEYWORDS:
;   WAVERANGE: Wavelength range that the ULYSS fitting procedure will use.
;              
;   VELSCALE:  Specifies the sampling step to be used by ULYSS, in km/s
;
;   /QUIET:    Verbosity control: suppress the messages
;
; SEE ALSO:
;     Definition of other components: ULY_SSP, ULY_TGM, ULY_LINE
;
; HISTORY:
;              Philippe Prugniel, 2008/05 created
;-
; CATEGORY:    ULY_UTIL
;------------------------------------------------------------------------------
function uly_star_init, cmp, WAVERANGE=lamrange, VELSCALE=velscale, QUIET=quiet

compile_opt idl2
on_error, 2

init_data = *cmp.init_data

if n_elements(*init_data.lsf_file) ne 0 then begin
    lsf = *init_data.lsf_file
endif

if size(init_data.file, /TYPE) eq 8 then begin
   file_in = uly_spect_logrebin(init_data.file, velscale, WAVERANGE=lamrange)
endif else begin
   if size(init_data.file, /TYPE) eq 7 then begin
      if file_test(init_data.file) ne 1 then begin
         message, 'Error, file does not exist', /INFO
         return, 0
      endif
      file_in = uly_spect_read(init_data.file, lamrange[0], lamrange[1], $
                               VELSCALE=velscale, QUIET=quiet)
   endif
endelse 

if uly_spect_get(file_in) ne 1 then begin
    if size(init_data.file, /TYPE) eq 7 then $
      message, 'Could not read the file ('+init_data.file+')', /INFO $
    else $
      message, 'Invalid input structure', /INFO 
    return, 0
endif


cmp.eval_fun = 'STAR'

;Convlove LOSVD in case giving lsf_file
if n_elements(lsf) gt 0 then begin
   uly_spect_lsfconvol, lsf, file_in
endif

cmp.eval_data = ptr_new(*file_in.data)

cmp.npix = (size(*file_in.data,/DIM))[0]
cmp.start = double(file_in.start)
cmp.step = double(file_in.step)
cmp.sampling = fix(file_in.sampling)
if n_elements(*file_in.goodpix) gt 0 then begin
    cmp.mask = ptr_new(bytarr(cmp.npix))
    (*cmp.mask)[*file_in.goodpix] = 1 
endif

uly_spect_free, file_in

return, cmp  

end

;----------------------------------------------------------------------------
;+
; NAME:        
;              ULY_STAR
;
; PURPOSE:     
;              Define a star component 
;
; USAGE:
;              cmp = uly_star(file_in, LSF=lsf_file, NAME=name, $
;              WL=lim_weight)
;
; ARGUMENT:
;   FILE_IN:   Input spectrum file name, spectrum structure (see 
;              ULY_SPECT_ALLOC) or array of file names or spectrum structures.
;              In the package we supply a solar spectrum as template, the users
;              can find and use their own template star(s).
;
; KEYWORDS:
;   LSF:       Name of the file containing a relative LSF to be injected 
;              in the template. 
;
;   NAME:      A string, name of this component. By default a unique name is  
;              generated.
;
;   WL:        Limits for the weight of each model component [min_w, max_w].
;              The weight is in data units over component units.
;              By default the weight is constrained to be positive, to
;              suppress any constraint, set: WL=[-1,1]*machar(/DOUBLE)).xmax
;              This constraint is ignored when ULySS fits a single component.
;
; DESCRIPTION:
;     Return a star fitting component or array of star component to be passed
;     to ULYSS. 
;     A star component has no free parameter. It can be a stellar template
;     or any sort of spectrum. It is typically used to fit the LOSVD.
;     Using an array of such components allows to search for the best positive
;     linear combination of spectra matchin an observation: this approach,
;     called 'optimal template matching', is commonly used to analyse
;     the kinematics of galaxies.
;
;     Note that a TGM component (see ULY_TGM) can be used to determine
;     in the same time the LOSVD and the atmospheric parameters that would
;     best represent the fitted spectrum.
;
; OUTPUT:     
;     Single star cmp struct, for more detailed explanation please
;     check ULY_FIT.
;
; REQUIRED FUNCTION:
;                 ULY_CMP_NAME
;
; EXAMPLE:
;              cmp = uly_star(uly_root+'/models/sun.fits')
;
; HISTORY:
;              Philippe Prugniel, 2008/05 created
;
;-
; CATEGORY:    ULY_STAR
;------------------------------------------------------------------------------
function uly_star, file_in,       $
                   LSF=lsf_file,  $
                   NAME=name,     $
                   WL=lim_weight

  compile_opt idl2
  on_error, 2

  if n_elements(file_in) gt 1 then $
     return, [uly_star(file_in[0],LSF=lsf_file,NAME=name,WL=lim_weight), $
              uly_star(file_in[1:*],LSF=lsf_file,NAME=name,WL=lim_weight)]

  if size(file_in, /TYPE) eq 0 then begin
     print, 'Usage: ULY_STAR, <filename|spect>, ...'
     mess = 'Error, no filename or spect was given'
     message, mess, /INFO
     return, 0
  endif else if size(file_in, /TYPE) eq 7 then begin
     if file_test(file_in) ne 1 then begin
        print, 'Usage: ULY_STAR, <filename|spect>, ...'
        mess = 'Error, file does not exist'
        if file_in ne '' then mess += ' ('+file_in+')'
        message, mess, /INFO
        return, 0
     endif
  endif else if size(file_in, /TYPE) eq 8 then begin
     if uly_spect_get(file_in) ne 1 then begin
        print, 'Usage: ULY_STAR, <filename|spect>, ...'
        message, 'Invalid input structure', /INFO
        return, 0
     endif
  endif
  namec = uly_cmp_name(name)
  s_wl = size(lim_weight)
  if s_wl[0] gt 0 then if not s_wl[0] eq 2 then $ 
     message, 'The weight limits have to be of the type arr(2)'

  lw = (n_elements(lim_weight) eq 0) ? (machar(/DOUBLE)).xmax*[0,1] : double(lim_weight)

  if size(file_in, /TYPE) eq 7 then descr = 'file:' + file_in + ' ' $
  else descr = 'spect'
  if n_elements(lsf_file) eq 1 then descr += 'lsf:' + lsf_file 

  init_data = {file:file_in, lsf_file:ptr_new(lsf_file)}

; return the structure describing this component
  return, {name:namec, $
           descr:descr, $
           init_fun:'uly_star_init', $
           init_data:ptr_new(init_data), $
           eval_fun:'', $
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
end


