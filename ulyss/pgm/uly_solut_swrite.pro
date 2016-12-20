;+
; NAME:
;              	  ULY_SOLUT_SWRITE
;
; PURPOSE:
;      Save the fit solution structure returned by ULY_FIT in a FITS file
;
; USAGE:
;		  uly_solut_swrite, <SOLUTION>, <FILE_OUT> [,STATUS=status]
;
; ARGUMENTS:
;   <SOLUTION>:  
;		  Structure created by ULY_FIT containing the fitting results.
;
;   <FILE_OUT>:  
;		  Name of the FITS file where the solution is saved.
;                                                            
; KEYWORD:
;   [STATUS=status]:  
;	          Returned status code. This code is either 0 for OK or an
;	          error number.
;
; DESCRIPTION:
;      A solut structure produced by ULY_FIT contains the solution of the
;      it as tabular data (handled by ULY_SOLUT_TREAD / ULY_SOLUT_TWRITE)
;      and as FITS data (original spectrum, best fit, multiplicative
;      polynomial...) handled by ULY_SOLUT_SREAD and the present routine.
;
; KNOWN BUGS:
;      This program cannot save 2d or 3d fits files
;
; HISTORY:
;                 Mina Koleva, 24/02/2006, created
;
;-
; CATEGORY:       ULY
;------------------------------------------------------------------------------

pro uly_solut_swrite, solution, file_out, STATUS=status

  status = 0

  if size(solution, /TYPE) ne 8 then begin 
     status = 1
     message, 'The solution is not a structure, can not write..', /INFO
     return
  endif

  if size(solution, /DIM) gt 1 then begin ; array of solut ... check consistency
     if min(solution.start) ne max(solution.start) then begin
        status = 1
        message, 'The solut array is not homogeneous (discrep. start)..', /INFO
        return
     endif
     if min(solution.step) ne max(solution.step) then begin
        status = 1
        message, 'The solut array is not homogeneous (discrep. step)..', /INFO
        return
     endif
  endif

  file_out2 = file_out+'.fits'  

; create output file
  fits_open, file_out2, fcb, /WRITE ; create a new FITS file
  h = uly_spect_get(solution[0], /HDR)
  if n_elements(h) le 1 then h = ['']

; add the spectral WCS keywords
  sxaddpar, h, 'CRVAL1', solution[0].start, 'ln(wave/0.1nm)'
  sxaddpar, h, 'CRPIX1', 1
  if solution[0].sampling eq 1 then begin
     sxaddpar, h, 'CDELT1', solution[0].step, 'ln(wave/0.1nm)'
     sxaddpar, h, 'CD1_1', solution[0].step, '0.1nm /in log(e)'
     sxaddpar, h, 'CTYPE1','AWAV-LOG','wavelength in air'
     sxaddpar, h, 'CUNIT1','Angstrom','unit for for the wavelength axis'
  endif else if solution[0].sampling eq 0 then begin
     sxaddpar, h, 'CDELT1', solution[0].step, '[0.1nm]'
     sxaddpar, h, 'CD1_1', solution[0].step, '[0.1nm]'
     sxaddpar, h, 'CTYPE1','AWAV','wavelength in air'
     sxaddpar, h, 'CUNIT1','Angstrom','unit for for the wavelength axis'
  endif

;add keywords to the header
  sxaddpar, h, 'FILENAME', solution[0].title,'name of the input file'
  sxaddpar, h, 'TITLE', solution[0].title,'name of the input file'
  sxaddpar, h, 'DOF_FACT', solution[0].dof_factor,'Degree of freedom factor'
  if tag_exist(solution[0], 'MODEL_ID') then $
     sxaddpar, h, 'MODEL_ID', solution[0].model_id,'model identification'

  sxdelpar, h, ['EXTNAME', 'BUNIT', 'DC-FLAG', 'WAT0_*', 'WAT1_*', 'WFITTYPE']

; A normal solut structure contains the following data arrays as members
; of the structure:
;   member     array name
;   DATA       DATA
;   ERROR      ERROR
;   BESTFIT    BESTFIT
;   MULCONT    MULTIPLICATIVE CONTINUUM
;   ADDCONT    ADDITIVE CONTINUUM
;   MASK       MASK
; However, any of them, except DATA, may be missing, but those which exist must
; have the same dimention, otherwise the structure is considered as invalid.
  
; I have to check all the tags, check the dimension of the data
; arrays, and store the data into 
  h1 = h
  sxaddpar, h1, 'FORMAT','solut', 'File written by ULYSS'
  
  i = 0
  if tag_exist(solution, 'data') then begin
     dim = size(*solution[0].data, /DIM)
     i += 1
     sxaddpar, h1, 'ARRAY'+strtrim(i,2),'DATA'
  endif
  if tag_exist(solution, 'err') then if n_elements(*solution.err) gt 0 then begin
     if i eq 0 then dim = size(*solution[0].err, /DIM) $
     else if array_equal(dim, size(*solution[0].err, /DIM)) eq 0 then status = 1
     i += 1
     sxaddpar, h1, 'ARRAY'+strtrim(i,2),'ERROR'
  endif
  if tag_exist(solution, 'bestfit') then begin
     if i eq 0 then dim = size(solution[0].bestfit, /DIM) $
     else if array_equal(dim, size(solution[0].bestfit, /DIM)) eq 0 then status = 1
     i += 1
     sxaddpar, h1, 'ARRAY'+strtrim(i,2),'BESTFIT'
  endif
  if tag_exist(solution, 'mulcont') then begin
     if i eq 0 then dim = size(solution[0].multcont, /DIM) $
     else if array_equal(dim, size(solution[0].mulcont, /DIM)) eq 0 then $
        status = 1
     i += 1
     sxaddpar, h1, 'ARRAY'+strtrim(i,2),'MULTIPLICATIVE CONTINUUM'
  endif
  if tag_exist(solution, 'addcont') then begin
     if i eq 0 then dim = size(solution[0].addcont, /DIM) $
     else if array_equal(dim, size(solution[0].addcont, /DIM)) eq 0 then $
        status = 1
     i += 1
     sxaddpar, h1, 'ARRAY'+strtrim(i,2),'ADDITIVE CONTINUUM'
  endif
  if tag_exist(solution, 'mask') then begin
     if i eq 0 then dim = size(solution[0].mask, /DIM) $
     else if array_equal(dim, size(solution[0].mask, /DIM)) eq 0 then $
        status = 1
     i += 1
     sxaddpar, h1, 'ARRAY'+strtrim(i,2),'MASK'
  endif

  if status eq 1 or n_elements(dim) eq 0 then begin
     status = 1
     message, 'The solut structure is inconsistent (arrays of different sizes)..', /INFO
     return
  endif
  iarray = i
  nsol = n_elements(solution)
  data = fltarr(product([dim,nsol]),iarray)

  i = 0
  if nsol eq 1 then begin
     if tag_exist(solution, 'data') then begin
        data[*,i] = (*solution[0].data) 
        i += 1
     endif
     if tag_exist(solution, 'err') then if n_elements(*solution.err) gt 0 then begin 
        data[*,i] = (*solution[0].err) 
        i += 1
     endif
  endif else begin
     if tag_exist(solution, 'data') then begin
        data = reform(data, product(dim), nsol, iarray, /OVER)
        for k=0,nsol-1 do data[*,k,i] = (*solution[k].data) 
        data = reform(data, product([dim,nsol]),iarray, /OVER)
        i += 1
     endif
     if tag_exist(solution, 'err') then if n_elements(*solution.err) gt 0 then begin 
        data = reform(data, product(dim),nsol, iarray, /OVER)
        for k=0,nsol-1 do data[*,k,i] = (*solution[k].err) 
        data = reform(data, product([dim,nsol]),iarray, /OVER)
        i += 1
     endif
  endelse   
  if tag_exist(solution, 'bestfit') then begin
     data[*,i] = solution.bestfit
     i += 1
  endif
  if tag_exist(solution, 'mulcont') then begin
     data[*,i] = solution.mulcont
     i += 1
  endif
  if tag_exist(solution, 'addcont') then begin
     data[*,i] = solution.addcont
     i += 1
  endif
  if tag_exist(solution, 'mask') then begin
     data[*,i] = solution.mask
     i += 1
  endif
  data = reform(data, [dim,nsol,iarray], /OVER)

;; write the fits file
  fits_write, fcb, data, h1, message=message
  fits_close, fcb

; 1) I think it would be more meaningful to write the various solout in
; different extensions : check in the procedure to fit long slit
; spectra what would appear to be the best.
; 2) We have to check how files with missing bestfit or other are handled  
; by uly_solut_sread, splot ...
; 3) I should prevent to write addcont and mulcont if they are not fitted
;   (maybe modify solut in order not to carry them)
; 4) We should study if it would not be better to store all the array
; as arrays rather than pointers for data and err
; 5) We may consider the possibility to have other ARRAY members ...

; cannot give EXTNAME to 1st HDU with fits_write ... update it after !
; (maybe we should not use fits_write)
;h = headfits(file_out2)         ; reread header
;sxaddpar, h, 'EXTNAME', 'FLUX'  ; Update EXTNAME value
;modfits, file_out2, 0, h, ext=0 ; Update extension hdr

end
