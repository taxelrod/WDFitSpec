;+
; NAME:                 ULY_SSP_WRITE
;
; PURPOSE:
;                       Write a SSP model grid to a FITS file
;
; USAGE:                uly_ssp_write, <model_grid>, <filename>
;
; ARGUMENTS:
;   <model_grid> :      Structure with the models (see ULY_SSP_READ)
;
;   <filename>   :      Name of the output FITS file
;
; DESCRIPTION:
;   Write the model structure <model_grid> into the file <filename> that
;   can be used later with ULY_SSP_READ or to define a component with ULY_SSP.
;
;   This function is in particular useful when a LSF has been injected in
;   an original model in order to use this new model to analyse
;   several spectra.
;
;   Description of the model file:
;   The data grid is stored as a N-dimensional array whose first dimension
;   is the wavelength. The other dimensions are the axes of the models:
;   AGE and METAL.
;   The coordinates on the wavelength axis are described with a standard
;   WCS, while for the other axese they are stored in extensions respectively
;   named AGE and METAL.
;   An additional extension, PARAM, if present, contains further information
;   on each spectrum of the model (like flux in given bands or Lick indices)
;   
;
; EXAMPLE:
;                       Read a grid, inject a LSF, write the modified grid
;     grid = uly_ssp_read(uly_root+'/models/PHR_Elodie31.fits', VEL=30.)
;     uly_spect_lsfconvol, 'lsf', grid
;     uly_ssp_write, grid, 'convolved_grid.fits'
;
;                       The LSF file 'lsf' used above can be generated as:
;     galaxy = uly_root+'/data/VazMiles_z-0.40t07.94.fits'
;     model = uly_ssp_extr(galaxy, $
;      uly_root+'/models/PHR_Elodie31.fits' ,[8000.,-0.4,0.], SIG=30)
;     cmp = uly_star(model)
;     uly_lsf, galaxy, cmp, 200, 100, FIL='lsf', /QUIET
;
; HISTORY:  
;           creation     Mina Koleva
;           update       Philippe Prugniel
;-
; CATEGORY:    ULY_SSP
;------------------------------------------------------------------------------
pro uly_ssp_write, model_grid, filename

; test the validity of arguments
if uly_ssp_chck (model_grid) ne 1 or size(filename, /TYPE) ne 7 then $
  print, 'Usage: ULY_SSP_WRITE, <model_grid>, <filename> '

if uly_ssp_chck (model_grid) ne 1 then message, '<model_grid> is not valid'
if size(filename, /TYPE) ne 7 then message, '<filename> is not valid'

; Initialize header for the output file
if n_elements(*model_grid.hdr) eq 0 then mkhdr, h,((*model_grid.data)[*,*,0]) $
else h = *model_grid.hdr

; Complete information in the Header
sxaddpar, h, 'HISTORY', 'model_grid created by ULY_SSP_WRITE'
sxaddpar, h, 'CRPIX1', 1, 'reference pixel'
sxaddpar, h, 'CRVAL1', model_grid.start, 'value at the center of ref pix'
sxaddpar, h, 'CD1_1', model_grid.step, 'step in wavelength'
sxaddpar, h, 'CDELT1', model_grid.step, 'step in wavelength'
if model_grid.sampling eq 0 then begin
    sxaddpar, h, 'CTYPE1','AWAV','wavelength in air, linear sampling'
endif else if model_grid.sampling eq 1 then begin
    sxaddpar, h, 'CTYPE1','AWAV-LOG','air wavelength, log sampling: ' $
      + string(model_grid.step*299798.458)+'km/s'
endif else message, 'do not handle yet sampling='+string(model_grid.sampling)

; New format 3/4D cube, plus extensions for ages and metallicities
; and eventually Mg/Fe.
; ages and metallicities are stored in two extensions: AGE and
; METALLICITY.
; If relevant, the Mg/Fe levels are stored in the keyword series MGFE
; Mask of goodpixels is stored in another extension: MASK

sxaddpar, h, 'CTYPE2', 'AGE', 'Axis 2 is age'
sxaddpar, h, 'CTYPE3', 'METAL','Axis 3 is metallicity'

if size(*model_grid.data, /N_DIM) eq 4 then begin
   for k=1,n_elements(*model_grid.o_mgfe) do $
      sxaddpar, h, 'MGFE'+strtrim(k,2), (*model_grid.o_mgfe)[k-1],'Mg/Fe on axis 4'
endif

sxdelpar, h, ['EXTNAME']

fits_open, filename, fcb, /WRITE ; create a new FITS file
fits_write, fcb, (*model_grid.data), h, MESSAGE=message
if message ne '' then message, message

sxdelpar, h, ['CTYPE1','CTYPE2','CTYPE3', $
              'CRVAL1','CRVAL2','CRVAL3', $
              'CRPIX1','CRPIX2','CRPIX3', $
              'CD1_1','CD2_2','CD3_3' $
             ]

fits_write, fcb, float(*model_grid.o_age), EXTNAME='AGE', h, MESSAGE=message
if message ne '' then message, message
fits_write, fcb, float(*model_grid.o_metal), EXTNAME='METAL', h, MESSAGE=message
if message ne '' then message, message

if n_elements(*model_grid.goodpix) gt 0 then begin
    sz0 = (size(*model_grid.data, /DIM))[0]
    msk = bytarr(sz0)
    msk[*model_grid.goodpix] = 1b
    fits_write, fcb, msk, EXTNAME='MASK', h, MESSAGE=message
    if message ne '' then message, message
endif

fits_close, fcb
; additional parameters, like fluxes in bands, or Lick indices may be
; stored in an attached bintable: ax1,ax2,ax3,...,para1,...

; A simple format is also to store a single bintable, one col beeing the
; spectra (stored each in one cell), and other cols age, met, ...
; It is clean and efficient, and many data specialists would prefer it
; But the grid structure (wcs) is not seen.

return
end
;--- end ---
