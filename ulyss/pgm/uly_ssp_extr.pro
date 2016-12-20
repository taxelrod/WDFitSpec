;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; NAME:                ULY_SSP_EXTR
;
; PURPOSE:             Interpolate in the model grid to extract a spectrum
;                      
; USAGE:               model = uly_ssp_extr(<galaxy>, <model_file>,
;                                           <pars>, SIGMA=sigma, /QUIET)
;
; DESCRIPTION:
;       Extract 1D spectrum from an interpolated grid of model spectra
;       at a given age, metallicity and [Mg/Fe]. Rebin it according to the  
;       WCS of the input spectrum <galaxy> and cut it in the wavelength range
;       of this latter spectrum.
;       
;       To simply interpolate a model at given age, metallicity and [Mg/Fe] 
;       first read the model grid with ULY_SSP_READ and after use the 
;       lower level function ULY_SSP_INTERP.
;     
;
; ARGUMENT:
;       <galaxy>      filename or spectrum structure serving as a model 
;                     for the WCS.
;       <model_file>  Identification of the model
;       <pars>        array of parameters of the model SSP ([age, met, mgf])
;                     
; KEYWORD:
;   SIGMA:   Gaussian velocity dispersion in km/s, to convolve the
;            model grid with when loading it.
;   QUIET:   Verbosity control.
;
; RETURNS:    
;       A 1D SSP spectrum as a spect structure, see ULY_SPECT_ALLOC
;            
; EXAMPLE:
;       To extract a model of PHR/Elodie with the same
;       SSP equvalent age/metallicity like Vazdekis/Miles's model 
;            model = $
;            uly_ssp_extr(uly_root+'/data/VazMiles_z-0.40t07.94.fits', $
;            uly_root+'/models/PHR_Elodie31.fits', [7980, -0.38, 0])
; 
; DEPENDENCE:
;            [uly_spect_read.pro], [uly_spect_alloc.pro],
;            [uly_spect_logrebin.pro],  uly_spect_get.pro,
;            uly_spect_free.pro, uly_ssp_read.pro, uly_ssp_interp.pro 
;
; HISTORY:
;            Creation  Mina Koleva 2008/03/11 
;-
; CATEGORY:  ULY_SSP
;-------------------------------------------------------------------------
function uly_ssp_extr, galaxy, model_file, pars, SIGMA=sigma, QUIET=quiet

; input checks
usage = 'Usage: spect = uly_ssp_extr(<galaxy>, <model_file>, <pars>)'
if n_elements(galaxy) eq 0 then begin
    print, 'Argument <galaxy> must be given'
    print, usage
    return, 0
endif
if size(model_file, /TYPE) ne 7 then begin
    print, 'Argument <model_file> must be a file name'
    print, usage
    return, 0
endif
if file_test(model_file) ne 1 then begin
    print, 'Argument <model_file> must be a file name'
    print, usage
    return, 0
endif
if n_elements(pars) eq 2 then pars = [pars,0] else $
if n_elements(pars) ne 3 then begin
    print, usage
    message, '"Pars" should be 2 or 3 elements array [age, [Fe/H], [Mg/Fe]]'
endif

if size(galaxy, /TYPE) ne 8 then begin
    spec = uly_spect_read(galaxy, QUIET=quiet)
endif else begin
    spec = uly_spect_alloc(SPECTRUM=galaxy)
endelse



;; define the lamrange in which to extract the model
lamrange =  uly_spect_get(spec, /WAVERANGE, STATUS=status)
if status ne 0 then begin
    message, 'Could not find WAVERANGE of spectrum', /INFO
    return, 0
endif

;; check if the observations are in log wl scale, if not -- rebin
if not (spec.sampling eq 1) then spec = uly_spect_logrebin(spec, /OVER)

velscale = spec.step * 299792.458d
if velscale eq 0 then message, 'Could not determine VELSCALE'

uly_spect_free, spec

pars[0]  = alog(pars[0])   ;age in log

;; read the grid of the models
model_grid = uly_ssp_read (model_file, $
                           WAVERANGE=lamrange, VELSCALE=velscale,        $
                           SIGMA=sigma, QUIET=quiet)

; -------------------------------------------------------------------------
;; evaluate the SSP  (this part is also in uly_evaluate)

if min(finite(pars[0])) eq 0 then message, 'Require a non-finite age'
if min(finite(pars[1])) eq 0 then message, 'Require a non-finite Fe/H'
if not size(*model_grid.data, /N_DIM) eq 4 then $   ;grid with Mg/Fe resolution
  pars[2] = dblarr(1) 

models = uly_ssp_interp(model_grid, double(pars))

; -------------------------------------------------------------------------
; Build the output spectrum structure

; Clean the input header
sxdelpar, *model_grid.hdr, ['NAXIS'+['2','3'], 'CTYPE'+['1','2','3'], 'CRPIX'+['1','2','3'], 'CRVAL'+['1','2','3'], 'CDELT'+['1','2','3'], 'CUNIT'+['1','2','3'], 'ULY_TYPE', 'CD1_1']

s = uly_spect_alloc(TITLE=model_grid.title, $
                    DATA=models, $
                    HEAD=*model_grid.hdr, $
                    START=model_grid.start, $
                    STEP=model_grid.step, $
                    SAMPLING=model_grid.sampling, $
                    GOODPIX=*model_grid.goodpix)

sxaddpar, *s.hdr, 'ULY_AGE', exp(pars[0]), '[Myr] Age extracted from the model grid'
sxaddpar, *s.hdr, 'ULY_FEH', pars[1], '[dex] Fe/H extracted from the model grid'

sxaddpar, *s.hdr, 'HISTORY', 'uly_ssp_extr '

uly_ssp_free, model_grid

return, s

end
