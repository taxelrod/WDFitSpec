;+
; NAME:    ULY_SSP_ALLOC ULY_SSP_FREE
;
; PURPOSE:
;          Allocate/free a structure containing the stellar population models
;
; USAGE:
;          grid = uly_ssp_alloc()
;
;          uly_ssp_free, grid
;
; DESCRIPTION:
;          ULY_SSP_ALLOC allocates an empty structure.
;          ULY_SSP_FREE frees an existing one.
;
;   This structure has the same tags as the spect structure 
;   (see ULY_SPECT_ALLOC) plus some additional ones. Therefore, many
;   functions designed to manipulate spect structures can be used, but
;   notice that (i) these functions will not handle the specific tags and
;   (ii) only ULY_SSP_FREE will correctly free the memory.
;
;   The anonymous structure contains the following tags:
;     title         string
;                Title, to be filled with the name of the grid                 
;     hdr        Array of character strings
;                FITS style header, to be filled with the hdr of the grid 
;     start      double precision
;                Wavelength of the centre of the first pixel
;                If 'sampling'=0, wavelength in Angstrom
;                If 'sampling'=1, log(wavelength)  (natural log)
;                If 'sampling'=2, wavelength can not be described with
;                start/step
;     step       double precision
;                Wavelength step between two pixels 
;                If 'sampling'=0, in Angstrom
;                If 'sampling'=1, logarithmic step
;                If 'sampling'=2, not a constant step
;     sampling   integer
;                Type of sampling scale
;                0 for linear in wavelength
;                1 for logarithmic scale  
;                2 for non-constant step
;     wavelen    Wavelength array, relevant if sampling=2
;     waverange  2 elements vector
;                The wavelength range for the data
;     goodpix    integer 1D array
;                List of valid pixels
;     data       real array of dimension 3 or 4
;                Pointer on Array of the values of the flux  
;     err        A pointer to null to make the structure compatible with 'spect'
;     o_age      Array of models' log(ages)
;     o_metal    Array of models' metallicities
;     o_mgfe     Array of models' [Mg/Fe], for a grid resolved in [Mg/Fe]
;     d2t        3 or 4D array containing the 2nd time derivatives
;                used for spline interpolation
;     dof_factor Dummy varaible to make the structure compatible with 'spect'
;
; RETURN:
;            empty model structure, to be filled with the (para)data of the
;            models 
;
; SEE ALSO:  ULY_SPECT_ALLOC
;
; AUTHOR:    Martin France, Mina Koleva, Philippe Prugniel
;            CRAL-Observatoire de Lyon
;
;-
; CATEGORY:  ULY_SSP
;------------------------------------------------------------------------------

; ***************************************************************************
; Deallocate a model structure
pro uly_ssp_free, model

compile_opt idl2
on_error, 0

if n_elements(model) ne 0 then $
  ptr_free, model.hdr, model.goodpix, model.data, model.err, $
  model.wavelen, model.waverange, $
  model.o_age, model.o_metal, model.o_mgfe, model.d2t
end

; ***************************************************************************
; Allocate a ssp model structure 
function uly_ssp_alloc

compile_opt idl2
on_error, 2

model_str = { title:'', $
              hdr:ptr_new(/ALLOCATE_HEAP),                     $
              start:1.0d, step:1.0d, sampling:1,              $
              wavelen:ptr_new(/ALLOCATE_HEAP),                 $
              waverange:ptr_new(/ALLOCATE_HEAP),               $
              goodpix:ptr_new(/ALLOCATE_HEAP),                 $
              data:ptr_new(/ALLOCATE_HEAP),                    $
              err:ptr_new(/ALLOCATE_HEAP),                     $
              o_age:ptr_new(/ALLOCATE_HEAP),                   $
              o_metal:ptr_new(/ALLOCATE_HEAP),                 $
              o_mgfe:ptr_new(/ALLOCATE_HEAP),                  $
              d2t:ptr_new(/ALLOCATE_HEAP),                     $
              dof_factor:1d                                    $
            }

return, model_str

end

;--- end ----------------------------------------------------------------------
