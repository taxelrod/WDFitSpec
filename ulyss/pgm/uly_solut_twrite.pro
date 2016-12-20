;+
; NAME:
;                 ULY_SOLUT_TWRITE
; PURPOSE:
;                 Save the solution of a fit to an ascii file
; USAGE:
;                 uly_solut_twrite, <solution>, <filename> [, HEAD='<command>']
;
; ARGUMENTS:
;   <solution>:   solut structure (or array) containing the solution(s).
;                 If <solution> is undefined, write only the header.
;   <filename>:   Name of the file where <solution> will be saved.
;                                                            
; KEYWORDS:
;    HEAD         Character string to write as header or comment. It may
;                 contain the command that generated <solution>.
;
; AUTHOR:
;                 Mina Koleva
; HISTORY:
;                 2008/02/05, MK, created.
;-
; CATEGORY:    ULY
;------------------------------------------------------------------------------
pro uly_solut_twrite, solution, file_res, HEAD=head

; NOTE:
;  We choosed to save the result in ASCII files because we did not
;  find a clean way to add records to existing FITS tables

usage = 'Usage: uly_solut_twrite, <solution>, <file_res>, HEAD=command'

if n_elements(file_res) eq 0 then begin
    message, 'No output file name...', /INFO
    print, usage
    return
endif
openw, lun, file_res, BUFSIZE=0, /GET_LUN, /APPEND
    
if n_elements(head) ne 0 then begin
;    printf, lun, '! Result file produced by the program ULYSS'
;    printf, lun, '! Date:' + systime()
;    printf, lun, '! Command: ' + head
;    printf, lun, '!
;    --------------------------------------------------------------'
    undefine, h
    if size(head, /TYPE) eq 7 then h = head
    get_date, dte, /TIMETAG
    sxaddpar, h, 'DATE', dte, 'When the header was written'
    h = ['BEGIN   ULYSS solut file (header section)',h]

    if n_elements(solution) gt 0 then begin
;      A change of the integer part of version number mark a full 
;      incompatibility (the reading pgm will not work or give wrong content)
;      A change in the decimal part indicate some incompatibility: in some
;      cases some information will not be read.
       sxaddpar, h, 'SOL_VERS', 2, 'ULySS solut file format version'
       if tag_exist(solution[0], 'data') then begin
          c = 299792.458d       ; [km/s] Speed of light
          wave_range = uly_spect_get(solution[0], /WAVERANGE, STATUS=status)
          if status eq 0 then begin
             sxaddpar, h, 'WAVE_MIN', wave_range[0], '[Angstrom] Start wavelength range'
             sxaddpar, h, 'WAVE_MAX', wave_range[1], '[Angstrom] End wavelength range'
          endif
       endif
       if tag_exist(solution[0], 'sampling') then $
          if solution[0].sampling eq 1 then $
             sxaddpar, h, 'VELSCALE', c*solution[0].step, '[km/s/pix] Velocity sampling'
       if tag_exist(solution[0], 'chisq') then $
          if solution[0].chisq eq 0 then $
             sxaddpar, h, 'NOERROR', 1, 'There was no error spectrum attached to DATA'
    endif

    printf, lun, '!'+h(where(strtrim(h) ne ''))
endif

if size(solution, /TYPE) ne 8 then begin
    if n_elements(solution) gt 0 then $
      message, 'The solution is not a structure, can not write..', /INFO
    close, lun
    free_lun, lun 
    return           ; nothing to write (or invalid)
endif

nsol = n_elements(solution)
cols_losvd = ['cz', 'e_cz', 's_cz', 'sig', 'e_sig', 's_sig', $
              'h3', 'e_h3', 's_h3', 'h4', 'e_h4', 's_h4', $
              'h5', 'e_h5', 's_h5', 'h6', 'e_h6', 's_h6']

for i = 0, nsol-1 do begin

;   Write the chisq col, even when not relevant, for compat with 
;   uly_solut_tread in ULySS v1.2 and earlier
    cols = 'chisq '
    vals = [0d]
    unts = ''

    if tag_exist(solution[i], 'snr') then begin
        cols += 'snr '
        vals = [vals, solution[i].snr]
    endif

    if tag_exist(solution[i], 'losvd') then begin
        moments = n_elements(solution[i].losvd)
        cols += strjoin(cols_losvd[0:3*moments-1] + ' ')

        for k = 0, moments-1 do $
          vals = [vals, solution[i].losvd[k], solution[i].e_losvd[k], solution[i].s_losvd[k]]
    endif

;   weight and params of the components
    if tag_exist(solution[i], 'cmp') then begin
        cmp = solution[i].cmp
        for j = 0, n_elements(cmp)-1 do begin   ; loop of each cmp

;           weight and e_weight of the component
            cols += ' ' + cmp[j].name + '_weight ' + $
              cmp[j].name + '_e_weight ' + cmp[j].name + '_l_weight'
            unts += ' -'   ; weight may have unit, but we do not store yet!
            vals = [vals, cmp[j].weight, cmp[j].e_weight, cmp[j].l_weight]

;           parameters and errors  
            n_para = -1
            if ptr_valid(cmp[j].para) then n_para = n_elements(*cmp[j].para) -1
            for k = 0, n_para do begin
                if strlen((*cmp[j].para)[k].name) gt 0 then begin
                    cols += ' ' + cmp[j].name + '_' + (*cmp[j].para)[k].name +$
                      ' ' + cmp[j].name + '_e_' + (*cmp[j].para)[k].name+ $
                      ' ' + cmp[j].name + '_s_' + (*cmp[j].para)[k].name+ $
                      ' ' + cmp[j].name + '_g_' + (*cmp[j].para)[k].name
                    vals = [vals, (*cmp[j].para)[k].value, $ 
                                  (*cmp[j].para)[k].error, $
                                  (*cmp[j].para)[k].status, $
                                  (*(*cmp[j].para)[k].guess)[0]]
                    if n_elements((*cmp[j].para)[k].unit) eq 0 then unts += ' -' $
                    else if strlen(strtrim((*cmp[j].para)[k].unit)) eq 0 then $
                      unts += ' -' $
                    else unts += ' ' + (*cmp[j].para)[k].unit
                    if n_elements((*cmp[j].para)[k].dispf) ne 0 then $
                      if strlen(strtrim((*cmp[j].para)[k].dispf)) gt 0 then $
                      unts += ','+(*cmp[j].para)[k].dispf
                endif                
            endfor

        endfor
    endif

    vals = vals[1:*]

    format = '(2a,/,2a,/,2a)'
    printf, lun, $
      '!cols:', cols, $
      '!unts:', unts, $
      string(solution[i].chisq), strjoin(string(vals)), FORMAT=format

endfor

close, lun
free_lun, lun

end
;--- end ----------------------------------------------------------------------
