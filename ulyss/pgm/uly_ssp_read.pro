;+
; NAME:           ULY_SSP_READ
;
; PURPOSE:        Read a grid of population model and sample it according to 
;                 specifications
;
; USAGE:          grid = uly_ssp_read( model_file,                   $
;                          TEMPLATE_GRID=template_grid,              $
;                          WAVERANGE=waverange, VELSCALE=velscale,   $
;                          AGERANGE=agerange, METALRANGE=metalrange, $
;                          LSF=lsf_file,                             $
;                          SIGMA=sigma, QUIET=quiet )
;
; ARGUMENTS:
;    model_file:  Name of the FITS file where the model are read
;
; INPUT KEYWORDS:
;    TEMPLATE_GRID (optional)
;        If given, ULY_SSP_READ checks that its value is a model
;        grid matching the description given through the other arguments
;        (see ULY_SSP_CHCK). If the grid match those specification
;        it is returned directly. Otherwise a new grid is read and
;        returned.
;
;     VELSCALE (optional)
;        Size of the pixels in km/s. For example VELSCALE=30.
;        If VELSCALE is not given, the models are not rebinned, ie.
;        they are as stored in the disk file (in linear or log scale).
;
;     WAVERANGE (optional)
;        Two element array giving the limits of the range to read, 
;        in Angstrom. For example WAVERANGE=[4000d,4786d].
;        By default the whole data range is read.
;
;     AGERANGE, METALRANGE (optional)
;        Extracted region of the grid. Ages in log(Myr) and metallicities
;        in dex.
;        By default the whole data range is read.
;        For example AGERANGE=alog([100,15000]), METALRANGE=[-1.4,+0.7]
;
;     LSF_FILE (optional)
;        Name of a file containing the LSF (see ULY_LSF and
;        ULY_SSP_LSFCONVOL) to be injected in the grid.
;
;     SIGMA (optional)
;        If given, the models are convolved by a Gaussian whose standard
;        deviation is SIGMA [km/s]
;
;     /QUIET (optional)
;        Verbosity control. 
;
; RETURN:
;   Structure containing the models grid. See ULY_SSP_ALLOC.
;   In case of error, return the scalar value -1.
;
; DESCRIPTION:
;        A grid of population model contains synthetic spectra at a
;        sample of ages, metallicity and possibly [Mg/Fe]. It is
;        used to generate a model at given age, metallicity and [Mg/Fe]
;        by interpolation.
;        The function ULY_SSP_READ returns a structure containing
;        the grid and all its description (wcs, sampling points in
;        age, metallicity ...).
;
; EXAMPLE: 
;   grid = uly_ssp_read(uly_root+'/models/PHR_Elodie31.fits', $
;                          WAVERANGE=[4000d, 5500d], VELSCALE=30., $
;                          AGERANGE=alog([12,20000]), METALRANGE=[-2.4,0.7])
;
;   grid = uly_ssp_read(uly_root+'/models/PHR_Elodie31.fits', VELS=30)
;
; AUTHORS:  Philippe Prugniel, Martin France, Mina Koleva
;
;-
; CATEGORY:    ULY_SSP
;------------------------------------------------------------------------------

function derive_3d, TMPARR=tmparr, GRID=grid
;   DERIVE_3D   (internal routine)
;   compute second derivatives of the interpolating function along the
;   2nd axis in the data cube (normally age). Use ULY_2DERIV 
compile_opt idl2, hidden

s = size(tmparr)

if (s[0] ne 3) then return, -1
if (n_elements(grid) ne s[2]) then return,-1

d2 = tmparr*0.0
for j=0,s[3]-1 do begin
    d2[*,*,j] = uly_2deriv(grid,tmparr[*,*,j])
endfor
; note: it is straightforward to suppress this for loop by making
; uly_2deriv work on the 3rd dimension. But it is not faster, probably
; because of the way the elements are accessed in memory

return, d2
end
;-----------------------------------------------------------------------------

function uly_ssp_read, model_file,$
                       TEMPLATE_GRID=template_grid, $
                       WAVERANGE=waverange, VELSCALE=velscale, $
                       AGERANGE=agerange, METALRANGE=metalrange, $
                       LSF=lsf_file, $
                       SIGMA=sigma, QUIET=quiet

compile_opt idl2
on_error, 0

if uly_ssp_chck(template_grid,                          $
                WAVERANGE=waverange, VELSCALE=velscale, $
                AGERANGE=agerange,                      $
                METALRANGE=metalrange,                  $
                MODEL_FILE=model_file, QUIET=quiet      $
               ) eq 1 then return, template_grid

bad = 0
if n_elements(model_file) eq 0 then begin
    str = 'MODEL_FILE, must be defined'
    bad = 1
endif
if n_elements(model_file) ne 1 then begin
    str = 'MODEL_FILE, must be a scalar or 1 element array'
    bad = 1
endif
if bad ne 1 then begin
    if file_test(model_file, /READ) ne 1 then begin
        str='MODEL_FILE, must be a file '
        if size(model_file, /TYPE) eq 7 then str += '('+ model_file+') '
        bad = 1
    endif                
endif

if bad eq 1 then begin
    print, 'ULY_SSP_READ: Usage:'
    print, '       grid = uly_ssp_read(model_file, $'
    print, '                 TEMPLATE_GRID=template_grid, $'
    print, '                 WAVERANGE=waverange, VELSCALE=velscale, $'
    print, '                 AGERANGE=agerange, METALRANGE=metalrange, $'
    print, '                 LSF=lsf_file, $'
    print, '                 SIGMA=sigma, QUIET=quiet)'
    if n_elements(str) gt 0 then message, /INFO, str
    return, -1
endif

if not keyword_set(quiet) then print, 'Read the model ', model_file
time0 = systime(1)


if not keyword_set(quiet) then print, 'read the file ', model_file
    
fits_open, model_file, fcb, /NO_ABORT, MESSAGE=messg
if messg ne '' then begin
    message, messg, /INFO
    return, -1
endif

if uly_ssp_chck(template_grid) eq 1 then grid = template_grid $
else grid = uly_ssp_alloc()

fits_read, fcb, *grid.data, *grid.hdr, MESSAGE=messg
if messg ne '' then begin
    message, messg, /INFO
    return, -1
endif
grid.start = double(string(sxpar(*grid.hdr, "CRVAL1")))
grid.step = double(string(sxpar(*grid.hdr, "CD1_1")))

if strtrim(sxpar(*grid.hdr, "CTYPE1")) eq 'AWAV' then grid.sampling = 0 $
else grid.sampling = 1

if n_elements(waverange) eq 2 then begin  ; check if WAVERANGE is valid
    wr = uly_spect_get(grid, /WAVERANGE)
    if wr[0] gt waverange[1] or wr[1] lt waverange[0] then begin
        message, 'Invalid WAVERANGE: '+strjoin(strtrim(waverange,2),','), /INFO
        return, -1
    endif
endif

if size(*grid.data, /N_DIM) eq 4 then begin ; grid with Mg/Fe resolution
    *grid.o_mgfe = sxpar(*grid.hdr, "MGFE*")
endif

fits_read, fcb, *grid.o_age, EXTNAME='AGE', MESSAGE=messg, /NO_ABORT
if messg ne '' then begin
    message, messg+'(AGE hdu)' , /INFO
    return, -1
endif

fits_read, fcb, *grid.o_metal, EXTNAME='METAL', MESSAGE=messg, /NO_ABORT
if messg ne '' then begin
    message, messg+'(METAL hdu)', /INFO
    return, -1 
endif

fits_read, fcb, mask, EXTNAME='MASK', MESSAGE=messg, /NO_ABORT
if messg ne '' then begin
    npix = (size(*grid.data, /DIM))[0]
    mask = bytarr(npix) + 1B
endif 

fits_close, fcb

masked = where(mask gt 0, cnt)
if cnt gt 0 then *grid.goodpix = masked $
else if not keyword_set(quiet) then $
  print, 'Ignore the MASK ... there is no valid data'

if n_elements(agerange) gt 0 then begin ; cut in age
    nage = where(*grid.o_age ge agerange[0] and *grid.o_age le agerange[1],cnt)
    if cnt gt 0 then begin
        *grid.data = (*grid.data)[*,nage,*]
        *grid.o_age = (*grid.o_age)[nage]
    endif else if not keyword_set(quiet) then $
      print, 'The demanded range:', strtrim(agerange,2), $
      'contain 0 elements of the grid range:',*grid.o_age, $
      '...read all the grid'
endif

if n_elements(metalrange) gt 0 then begin ; cut in metallicity
    nmet = where(*grid.o_metal ge metalrange[0] and *grid.o_metal le metalrange[1], cnt)
    if cnt gt 0 then begin
        *grid.data = (*grid.data)[*,*,nmet]
        *grid.o_metal = (*grid.o_metal)[nmet]
    endif else if not keyword_set(quiet) then $
     print, 'The demanded range:', strtrim(metalrange,2), $
      'contain 0 elements of the grid range:',*grid.o_metal, $
      '...read all the grid'
endif

grid.title = model_file
if n_elements(waverange) ne 0 then *grid.waverange = waverange

log_rebin = 0
if n_elements(velscale) eq 1 or n_elements(lsf_file) eq 1 or $
   n_elements(sigma) eq 1 then log_rebin = 1

if log_rebin then begin ;   rebin in log wavelength
    grid = uly_spect_logrebin(grid, velscale, WAVERANGE=waverange, /EXACT, /OVERWRITE)
    if not uly_spect_get(grid, /VALID) then begin
        message, 'Could not log-sample the grid', /INFO
        return, -1
    endif
endif

if not keyword_set(quiet) then print, 'model reading time=', systime(1)-time0

if n_elements(sigma) gt 0 then begin
    if grid.sampling ne 1 then message, $
       'ULY_SSP_READ: cannot convolve because the grid is not log-sampled'
    if not keyword_set(quiet) then begin
        print,'Convolve the model with Gaussian, sigma=', sigma, ' km/s'
    endif
    uly_slit, 0d, sigma, 0d, 0d, velscale, lsf
    *grid.data = convol((*grid.data), lsf, /EDGE_TRUNCATE)
endif

if n_elements(lsf_file) ne 0 then begin
    uly_spect_lsfconvol, lsf_file, grid
endif

timed0 = systime(1)

if size(*grid.data, /N_DIM) le 3 then begin ; normal grid
    *grid.d2t = derive_3d(TMPARR=*grid.data, GRID=*grid.o_age)
endif else begin                ; grid with Mg/Fe resolution
    form = size(*grid.data,/DIM)
    tot = n_elements((*grid.data)[*,*,*,0])
    *grid.d2t = [reform(derive_3d(TMPARR=(*grid.data)[*,*,*,0], GRID=*grid.o_age), tot),$
                 reform(derive_3d(TMPARR=(*grid.data)[*,*,*,1], GRID=*grid.o_age),tot)]
    *grid.d2t = reform(*grid.d2t, form)
endelse

if not keyword_set(quiet) then print, 'model deriving time=', systime(1)-timed0

return, grid

end
