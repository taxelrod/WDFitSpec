;+
; NAME:           ULY_SSP_INIT
;
; PURPOSE:        Initialize a SSP component
;
; USAGE:          cmp_new = uly_ssp_init(cmp, 
;                                        WAVERANGE=waverange,
;                                        VELSCALE=velscale, 
;                                        QUIET=quiet)
;
; DESCRIPTION:    This initialization function is automatically
;                 executed by ULY_FIT_INIT when a component has been
;                 defined using ULY_SSP (fit of a Single Stellar Population).
;
;                 This function is recorded in the member:
;                 cmp.init_fun of the component structure.
;
;                 The initialization consist in reading the grid of
;                 SSPs (ULY_SSP_READ) sample it according to WAVERANGE
;                 and VELSCALE and possibly inject the LSF.
; ARGUMENTS:
;    cmp:         a component defined using ULY_SSP (look its doc for
;                 more details)
;
; KEYWORDS:
;    WAVERANGE:   wavelength range for which the models are read
;    VELSCALE:    if this keyword is set the models are going to be
;                 rebined in log with a step in km/s equal to the velscale
;                 value
;    QUIET:       verbosity control 
;
; DEPENDENCE:     ULY_SSP_READ.pro
;
; HISTORY:        creation PhP
;-
; CATEGORY:       ULY_SSP
;----------------------------------------------------------------------------
; initialize a fit component for a SSP
function uly_ssp_init, cmp, WAVERANGE=lamrange, VELSCALE=velscale, QUIET=quiet

compile_opt idl2
on_error, 2

init_data = *cmp.init_data

template_grid = uly_ssp_read(init_data.model,$
                             TEMPLATE_GRID=*init_data.template_grid,     $
                             WAVERANGE=lamrange, VELSCALE=velscale,      $
                             LSF=*init_data.lsf_file,                    $
                             QUIET=quiet)

if size(template_grid, /TYPE) ne 8 then begin
    str=''
    if size(init_data.model, /TYPE) eq 7 then str = '('+ init_data.model[0]+') '
    message, 'Could not read SSP '+str, /INFO
    return, cmp
endif

*init_data.template_grid = template_grid

s1 = size(*template_grid.data)
s4 = size(*template_grid.d2t)
if (s1[0] gt 4) then begin
    message, 'Model grid has invalid dimensions', /INFO
    cmp.sampling = -1
    return, cmp                    ; an error occured
endif
if not array_equal(s1,s4) then begin
    message,$
      'The 2nd deriv. of age and templates must have the same size/type' , /INF
    cmp.sampling = -1
    return, cmp                    ; an error occured
endif

; load the parameters and data in the comp struct

agemin = min(*template_grid.o_age,max=agemax) ; log(age)
metmin = min(*template_grid.o_metal,max=metmax)
if size(*template_grid.data, /N_DIM) eq 4 then begin ;  grid with Mg/Fe resol
   inpmgfe = *template_grid.o_mgfe 
   (*cmp.para)[2].name = 'Mg/Fe'
   (*cmp.para)[2].unit = 'dex'
endif else inpmgfe = [0.,0.]
mgfemin=min(inpmgfe,max=mgfemax)

if total((*cmp.para)[0].limits ne [0d,0d]) eq 0 then $
   (*cmp.para)[0].limits = [agemin,agemax] $ ; log age limits
else if $
    (*cmp.para)[0].limits[0] lt agemin or (*cmp.para)[0].limits[1] gt agemax $
then begin
    message, /INFO, $
       'The required AGE limits (' + string([agemin,agemax],FOR='(a,1h,,a)')+ $
      ') are out of the model bounds (' + $
      string((*cmp.para)[0].limits,FOR='(a,1h,,a)') + ')'
    cmp.sampling = -1
    return, cmp                    ; an error occured
endif

if total((*cmp.para)[1].limits ne [0d,0d]) eq 0 then $
  (*cmp.para)[1].limits = [metmin, metmax] $ 
else if $
    (*cmp.para)[1].limits[0] lt metmin or (*cmp.para)[1].limits[1] gt metmax $
then begin
    message, /INFO, $
       'The required MET limits ('+strtrim((*cmp.para)[1].limits[0],2)+','+strtrim((*cmp.para)[1].limits[1],2)+') are out of the model bounds ('+strtrim(metmin,2)+','+strtrim(metmax,2)+')' 
    cmp.sampling = -1
    return, cmp                    ; an error occured
endif

if total((*cmp.para)[2].limits ne [0d,0d]) eq 0 then $
  (*cmp.para)[2].limits = [mgfemin-0.1, mgfemax+0.3]

cmp.eval_fun = 'SSP'
if ptr_valid(cmp.eval_data) then *cmp.eval_data = template_grid $
else cmp.eval_data = ptr_new(template_grid)

cmp.npix = (size(*template_grid.data,/DIM))[0]
cmp.start = double(template_grid.start)
cmp.step = double(template_grid.step)
cmp.sampling = fix(template_grid.sampling)
if n_elements(*template_grid.goodpix) gt 0 then begin
   if ptr_valid(cmp.mask) then *cmp.mask = bytarr(cmp.npix) $
   else cmp.mask = ptr_new(bytarr(cmp.npix))
   (*cmp.mask)[*template_grid.goodpix] = 1 
endif

return, cmp  ; OK, no error

end

;----------------------------------------------------------------------------
;+
; NAME:           ULY_SSP
;
; PURPOSE:        Define a SSP component 
;              
; USAGE:
;                 cmp = uly_ssp(MODEL_FILE=model_file, 
;                               DATA=template_grid,  
;                               LSF=lsf_file,                           
;                               AL=a_limits, 
;                               ZL=z_limits, 
;                               MGFEL=mg_limits,
;                               AG=a_guess, 
;                               ZG=z_guess, 
;                               MGFEG=mg_guess, 
;                               NAME=name, 
;                               FIXPAR=fixpar, 
;                               WL=lim_weight)
;
;
; DESCRIPTION:    
;      Return a CMP structure describing a SSP component that can be passed to 
;      the program ULYSS.
;
;      A SSP is a function of age, metallicity and possibly of [Mg/Fe]. The set of 
;      models is defined by MODEL_FILE, by default Pegase-HR models computed with 
;      ELODIE.3.1 and Salpeter IMF will be used.
;      A list of the models included in the basic distribution of the package
;      can be found on the models page. And other models can be downloaded
;      from the same page.
;
;      To define a n-bursts model use n times ULY_SSP and make an array
;      of these cmp_n: cmp = [cmp_1, cmp_2,... cmp_n]
;
;      Note that the SSP components are not obligarory SSPs, they may for example 
;      be population with constant SFR. They just depend at maximum on age, metallicity 
;      and [Mg/Fe]. (An exponentially decreasing SFR cannot be modeled with this CMP.
;
; KEYWORDS:
;   MODEL_FILE:  Name of the FITS file containing the SSP grid.
;                Defaulted to models/PHR_Elodie31.fits (PEGASE_HR
;                grid with Elodie.3.1 library).
;                An alternative choice, provided in the distribution
;                is to use the Vazdekis models build with the Miles
;                library (uly_root+'/models/Vaz_Miles.fits').
;
;   DATA:        This keyword is used for template model (input and output). 
;                Normally ULYSS reads the model grid specified with MODEL_FILE.
;                If <template_grid> is the correct model for the required 
;                action, then it uses it without reloading. 
;                The same keyword is used for output of the loaded model, ie.,
;                you can set this keyword to an empty named variable where 
;                the loaded model will be stored.                
;
;   LSF:         Name of the file containing a relative LSF to be injected in 
;                the template. 
;
;   AL:          [Myr] default=limits of the grid
;                Limits for the age. dblarr(2). 
;
;   ZL:          [dex] default=limits of the grid
;                dblarr(2). Limits for the metallicities.
;
;   MGFEL:       [dex] default=limits of the grid
;                dblarr(2); Limits for [Mg/Fe].
;
;   AG:          [Myr] default=8000
;                Guess for the age. This can be a scalar or an array.
;
;   ZG:          [dex] default=-0.4
;                Guess for the metallicity. This can be a scalar or an array.
;
;   MGFEG:       [dex] default=0.0
;                Guess for the [Mg/Fe]. This may be an array. 
;
;   NAME:        Character string to name the component. By default a unique
;                name is generated.
;
;   FIXPAR:      array[3], for [age, Fe/H, Mg/Fe]
;                the elements which have to be fixed 0/1 - free/fixed
;
;   WL:          limits for the weight of each model component [min_w, max_w]
;                The weight is in data units over component units.
;                By default the weight is constrained to be positive, to
;                suppress any constraint, set: WL=[-1,1]*machar(/DOUBLE)).xmax
;                This constraint is ignored when ULySS fits a single component.
;
; EXAMPLE:
;                To fit Vaz-Miles SSP with PHR-ELODIE
;                 cmp=uly_ssp()
;                 ulyss, uly_root+'/data/VazMiles_z-0.40t07.94.fits', cmp
;
; HISTORY:
;     CREATION    2008/07/24 Mina Koleva
;-
; CATEGORY:       ULY_SSP
;------------------------------------------------------------------------------

; define a fit component for a SSP
FUNCTION uly_ssp, MODEL_FILE=model_file, DATA=template_grid,     $
                  LSF=lsf_file,                              $
                  AL=a_limits, ZL=z_limits, MGFEL=mg_limits, $
                  AG=a_guess, ZG=z_guess, MGFEG=mg_guess,    $
                  NAME=name, FIXPAR=fixpar, WL=lim_weight

compile_opt idl2

on_error, 2

common uly_path, uly_root
if n_elements(model_file) eq 0 then begin
    if n_elements(template_grid) gt 0 then modelf = template_grid.title $
    else modelf = uly_root+'/models/PHR_Elodie31.fits'
endif else modelf = model_file

if n_elements(a_guess) eq 0  then a_guess = 8000
if n_elements(z_guess) eq 0  then z_guess = -0.4
if n_elements(mgfe_guess) eq 0 then mgfe_guess = 0.0

s5 = size(a_limits)
if s5[0] gt 0 then $
  if not s5[0] eq 2 then $
  message, 'The age limits have to be of the type arr(2)'

s6 = size(z_limits)
if s6[0] gt 0 then $
  if not s6[0] eq 2 then $
  message, 'The metallicity limits have to be of the type arr(2)'

s7 = size(lim_weight)
if s7[0] gt 0 then $
  if not s7[0] eq 2 then $ 
  message, 'The weight limits have to be of the type arr(2)'

; create the structure describing this component
init_data = {model:modelf, template_grid:ptr_new(template_grid), $
             lsf_file:ptr_new(lsf_file)}
descr = ''
if n_elements(modelf) eq 1 then descr += 'model:' + modelf + ' ' 
if n_elements(lsf_file) eq 1 then descr += 'lsf:' + lsf_file 

namec = uly_cmp_name(name)

cmp = {name:namec, $
       descr:descr, $
       init_fun:'uly_ssp_init', $
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
       lim_weig:(machar(/DOUBLE)).xmax * [0,1] $
      }

; load the parameters and data in the comp struct

*cmp.para = replicate({name:'', unit:'', guess:ptr_new(), step:1D-2, $
                       limits:[0d,0d], limited:[1,1], fixed:0S, $
                       value:0D, error:0D, dispf:''}, $
                      3)

(*cmp.para)[0].name = 'age'
(*cmp.para)[0].unit = 'Myr'
(*cmp.para)[1].name = 'Fe/H'
(*cmp.para)[1].unit = 'dex'

(*cmp.para)[0].dispf = 'exp'  ; says to display exp(para)

if n_elements(a_limits) gt 0 then (*cmp.para)[0].limits = alog(a_limits)
if n_elements(z_limits) gt 0 then (*cmp.para)[1].limits = z_limits
if n_elements(mg_limits) gt 0 then (*cmp.para)[2].limits = mg_limits

for k=0,2 do (*cmp.para)[k].guess = ptr_new((*cmp.para)[k].limits[0])
if n_elements(a_guess) ge 1 then $
  *(*cmp.para)[0].guess = alog(double(string(a_guess)))
if n_elements(z_guess) ge 1 then $
  *(*cmp.para)[1].guess = double(string(z_guess))
if n_elements(mg_guess) ge 1 then $
  *(*cmp.para)[2].guess = double(string(mg_guess))


(*cmp.para).step = 0.01d  ; step in log age;  ~10Myr on 1Gyr; O.O1 for all

(*cmp.para)[2].fixed = 1s ; fixed by default
if n_elements(fixpar) gt 0 then $
  (*cmp.para)[0:n_elements(fixpar)-1].fixed = fixpar

if n_elements(lim_weight) gt 0 then $
  cmp.lim_weig = double(lim_weight)

return, cmp

end

;-- end ----------------------------------------------------------------------

