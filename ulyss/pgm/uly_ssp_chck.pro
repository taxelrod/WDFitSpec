;+
; NAME:
;         ULY_SSP_CHCK
;
; PURPOSE:
;         Check if a model grid matchs specifications given in keywords
;
; USAGE:  
;         status = uly_ssp_chck(model_grid,                    $
;                     WAVERANGE=waverange, VELSCALE=velscale,  $
;                     AGERANGE=agerange, METALRANGE=metalrange,$
;                     MODEL_FILE=model_file, QUIET=quiet       $
;                    
;
; DESCRIPTION:
;         Check if a grid of population models previously read
;         (likely with ULY_SSP_READ) corresponds to some reading
;         specification. This may save to reload it.
;         This routine check if the model (MODEL_FILE, for example
;         'PHR_Elodie31.fits' or 'Vaz_Miles.fits'), WCS
;         (waverange, velscale) or range of parameters
;         (agerange, metalrange) match.
;         This routine compares the characteristics of the model grid
;         as they are stored in the structure (model_grid) with those
;         specified in keywords
;
;         If a model grid is passed in entry to ULY_SSP_READ through
;         the keyword MODEL_FILE, the present routine is called
;         to test if the grid is consistent with the other parameters.
;
; ARGUMENT:
;    model_grid: input  (not modified by this routine)
;         Structure created by eg. uly_ssp_read and containing the
;         model.
;         If not give, ULY_SSP_CHCK returns 0.
;
; KEYWORD:
;    The keywords are similar to those of ULY_SSP_READ
;     WAVERANGE: (optional)
;        Two element array giving the limits of the range to read, 
;        in Angstrom. For example WAVERANGE=[4000d,4786d].
;        By default the whole data range is read.
;     VELSCALE: (optional)
;        Size of the pixels in km/s. For example VELSCALE=30.
;        If VELSCALE is not given, the models are not rebinned, ie.
;        they are as stored in the disk file (in linear or log scale).
;     AGERANGE, METALRANGE: (optional)
;        Extracted region of the grid. Ages are in log(Myr) and metallicities
;        in dex.
;        By default the whole data range is read.
;        For example AGERANGE=alog([100,15000]), METALRANGE=[-1.4,+0.7]
;     MODEL_FILE: (string) Name of the file containing the grid of SSPs.
;     QUIET:  verbosity control
;
; RETURN:
;    integer (model_reload): 1 if (model_grid) does not corresponds to
;                              the specifications
;                            0 if (model_grid) corresponds to the specs
;
; EXAMPLE:
;        Read a grid of Pegase-HR/Elodie with given wavelength, age
;        and metallicity ranges:
;             grid = uly_ssp_read(uly_root+'/models/PHR_Elodie31.fits', $
;                          WAVERANGE=[4000d, 5500d], VELSCALE=30., $
;                          AGERANGE=alog([20,2000]), METALRANGE=[-1.4,-0.1])
;        Check it if <grid> corresponds to the reading specifications:
;             if uly_ssp_chck(grid,                                 $
;                     WAVERANGE=[4000d,5500d], VELSCALE=30,         $
;                     AGERANGE=[20, 2000], METALRANGE=[-1.4,-0.1],  $
;                     MODEL_FILE='models/PHR_elodie31.fits', /QUIET $
;                    ) eq 1 then print, 'It is the correct grid!'
;
; DEPENDENCE:
;            ULY_SSP_ALLOC.PRO
;
; AUTHOR:  Antoine Bouchard, Philippe Prugniel
;
;-
; CATEGORY:    ULY_SSP
;-----------------------------------------------------------------------------
function uly_ssp_chck, model_grid, $
                       WAVERANGE=waverange, VELSCALE=velscale, $
                       AGERANGE=agerange, METALRANGE=metalrange, $
                       MODEL_FILE=model_file,$
                       QUIET=quiet

compile_opt idl2
on_error, 2

if n_elements(model_grid) eq 0 then begin
    return, 0
endif

; Is it a structure?
if (size(model_grid, /structure)).type ne 8 then begin
  if not keyword_set(quiet) then $
    print, 'MODEL_GRID is not a valid structure: reloading grid'
  return, 0
endif

; Is this structure of the right type?
a = uly_ssp_alloc()
if total(tag_names(model_grid) ne tag_names(a)) gt 0 then begin
    uly_ssp_free, a
    if not keyword_set(quiet) then $
      print, 'MODEL_GRID is not a valid template: reloading grid'
    return, 0
endif
uly_ssp_free, a
    
if model_grid.sampling ne 1 and n_elements(velscale) then begin
    if not keyword_set(quiet) then $
      print, 'The loaded template is not in log scale: reloading grid'
    return, 0
endif

if n_elements(model_file) ne 0 then if model_grid.title ne model_file then begin
    if not keyword_set(quiet) then $
      print, 'The specified MODEL_FILE does not match with the one of ', $
      'the loaded template: reloading grid: ',model_grid.title, ' vs. ', model_file
    return, 0
endif 

if n_elements(waverange) ne 0 then $
  if total(*model_grid.waverange ne waverange) gt 0 then begin
    if not keyword_set(quiet) then $ 
      print, 'The wavelength range of the data does not match with ', $
      'the loaded template: reloading grid' 
    return, 0
endif 

if n_elements(velscale) ne 0 then $
  if 299792.458d*model_grid.step ne velscale then begin
    if not keyword_set(quiet) then $ 
      print, 'The specified velocity scale does not match with that of ', $  
      'the loaded template: reloading grid'
    return, 0
endif 

return, 1   ; the grid match the specs

end
