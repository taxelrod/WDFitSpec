;+
; NAME:
;                 ULY_SOLUT_SREAD
;
; PURPOSE:
;                 Read a solution file written by ULY_SOLUT_SWRITE
;
; USAGE:
;	          SOLUTION = uly_solut_sread(<solut_file>[, STATUS=<status>])
;
; ARGUMENT:
;   <solut_file>: Name of the FITS file containing the solution structure.
;                                                            
; KEYWORD:
;   [STATUS=<status>]:  
;	          Returned status code. This code is either 0 for OK or an
;	          error number.
;                                                            
; OUTPUT:
;   SOLUTION:
;	          Solution structure containing the results of a previous fit
;	          that has been saved on disk in the form of a FITS file.
;
; DESCRIPTION:
;   The solution is the output from ULY_FIT, and it written by 
;   ULY_SOLUT_SWRITE.
;
;   The solut FITS file is in the SDSS style (see ULY_SPECT_READ). It
;   consists in a image-type array of dimension 2. Thefirst axis is
;   the wavelength described by a standard WCS. The second axis is described
;   by the ARRAYi keyword series (ARRAY1 describe the first scan ...):
;           ARRAY1             'DATA    ' 
;           ARRAY2             'ERROR   ' 
;           ARRAY3             'BESTFIT ' 
;           ARRAY4   'MULTIPLICATIVE CONTINUUM' 
;           ARRAY5   'ADDITIVE CONTINUUM' 
;           ARRAY6             'MASK    ' 1:good, 0:bad
;  
;   The format of the output FITS file must be described here.
;
; AUTHORS:
;                 Antoine Bouchard, Philippe Prugniel
;
; HISTORY:
;                 2008/02/05
;-
; CATEGORY:       ULY
;------------------------------------------------------------------------------
function uly_solut_sread, solut_file, STATUS=status

status = 0

if (size(solut_file))[0] ne 0 or (size(solut_file))[1] ne 7 then begin
    status = 1
    message, 'Invalid file name', /INFO
    return, 0
endif

filen = solut_file+".fits"
if file_test(filen) eq 0 then filen = solut_file
if file_test(filen) eq 0 then begin
    status = 1
    message, 'File not found', /INFO
    return, 0
endif

fits_open, filen, fcb, /NO_ABORT, MESSAGE=messg
if messg ne '' then begin
    status = 1
    message, messg, /INFO
    return, 0
endif

fits_read, fcb, data, hdr, MESSAGE=messg
fits_close, fcb
if messg ne '' then begin
    status = 1
    message, messg, /INFO
    return, 0
endif
if sxpar(hdr, "NAXIS") eq 0 then begin
    message, 'data extension is empty', /INFO
endif

format = strtrim(sxpar(hdr, "FORMAT"))
if format ne 'solut' then message, 'This is not a solut file'

array = strtrim(sxpar(hdr, "ARRAY*"))

ndim = size(data,/N_DIM)
dim = size(data,/DIM)
narray = n_elements (array)
if dim[ndim-1] ne narray then begin
    message, /INFO, $
      'invalid file, ARRAY* keywords should describe the last axis of data '
endif
dim = dim[0:ndim-2]

tot = 1
for i=0, ndim-2 do tot *= dim[i]
data = reform(data, tot, narray, /OVER)

crv    = double(string(sxpar(hdr, "CRVAL1")))
cdv    = double(string(sxpar(hdr, "CD1_1")))
crp    = double(string(sxpar(hdr, "CRPIX1")))

ctype = sxpar(hdr,'CTYPE1')
if strmid(ctype,5,3) eq 'LOG' then sampling = 1 else sampling = 0

title      = sxpar(hdr, "TITLE")
dof_factor = sxpar(hdr, "DOF_FACT")
model_id   = sxpar(hdr, "MODEL_ID")

for i=0,n_elements(array)-1 do begin
    case array[i] of
        'DATA': flux = reform(data[*,i], dim)
        'ERROR': err = reform(data[*,i], dim)
        'BESTFIT': bestfit = reform(data[*,i], dim)
        'MULTIPLICATIVE CONTINUUM': mulcont = reform(data[*,i], dim)
        'ADDITIVE CONTINUUM': addcont = reform(data[*,i], dim)
        'MASK': mask = byte(reform(data[*,i], dim))
    endcase
endfor

return, { $
          title:title,                    $
          hdr:ptr_new(hdr),               $
          start:crv-(crp-1d)*cdv,         $
          step:cdv,                       $
          sampling:sampling,              $
	  refpix:crp, 		          $
          data:ptr_new(flux),             $
          err:ptr_new(err),               $
          bestfit:bestfit,                $
          mask:mask,                      $
          dof_factor:dof_factor,          $
          addcont:addcont,                $
          mulcont:mulcont,                $
          model_id:model_id               $              
}

end
