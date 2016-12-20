;+
; NAME:
;                 ULY_SOLUT_TPRINT
;
; PURPOSE:
;                 Print the solution of a fit
;
; USAGE:      
;                 uly_solut_tprint, solution     
;
; ARGUMENT:
;   solution:     Solut structure (or array of solut) or name of a file
;                 containing the solution(s)
;                                                            
; OUTPUT:
;                 Solut structure
;
; AUTHOR:
;                 Mina Koleva
;
; HISTORY:
;                 2008/02/05
;
;-
; CATEGORY:       ULY
;------------------------------------------------------------------------------
pro uly_solut_tprint, solution

if size(solution, /TYPE) eq 7 then begin
    solution = uly_solut_tread(solution) 
    readsol = 1
endif else readsol = 0

if size(solution, /TYPE) ne 8 then begin
  message, 'The solution is not a structure, cannot print..', /CONTINUE
  return
endif

for i = 0, (size(solution))[1]-1 do begin

   if tag_exist(solution[i], 'data') and $
      tag_exist(solution[i], 'mask') then $
         print, "Number of pixels used for the fit", $
                n_elements(where(solution[i].mask eq 1))

   if tag_exist(solution[0], 'losvd') eq 1 and $
      tag_exist(solution[0], 'e_losvd') then begin
      losnm = ["cz             :", $
               "dispersion     :", $
               "h3             :", $
               "h4             :", $
               "h5             :", $
               "h6             :"]
      losun = [" km/s", " km/s", " ", " ", " ", " "]
      print, ""
      for k=0,(size(solution[i].losvd, /DIM))[0]-1 do begin
         if solution[i].s_losvd[k] eq 0 then                            $
            print, losnm[k], solution[i].losvd[k], " +/- ",             $
                   strtrim(solution[i].e_losvd[k],2), losun[k]          $
         else if solution[i].s_losvd[k] eq 1 then                       $ 
            print, losnm[k], solution[i].losvd[k], losun[k], ' (fixed)' $
         else if solution[i].s_losvd[k] eq 2 then                       $
            print, losnm[k], solution[i].losvd[k],                      $
                   losun[k],' (pegged at low bound)'                    $
         else if solution[i].s_losvd[k] eq 3 then                       $
            print, losnm[k], solution[i].losvd[k],                      $
                   losun[k],' (pegged at high bound)'
      endfor
   endif
    
    print,'-----------------------------------------------'
    if solution[i].chisq gt 0 then $
       print, "chi^2          :", solution[i].chisq
    if solution[i].snr gt 0 then $
       print, "estimated SNR  :", solution[i].snr
    print,'-----------------------------------------------'
    
    if tag_exist(solution[i], 'cmp') then begin

        cmp = solution[i].cmp

        tw = total(cmp.l_weight) / 100
        
        nocontrib = ''
        for j=0,n_elements(cmp)-1 do begin
            if cmp[j].weight eq 0 then begin
               nocontrib += strtrim(j,2) + '(' + cmp[j].name + ') '
               continue
            endif
            print, 'cmp #', strtrim(j,2), '  ', cmp[j].name
            if n_elements(cmp) gt 1 then begin
               print, 'Light fraction :', cmp[j].l_weight / tw, ' +/- ', $
                 cmp[j].e_weight / cmp[j].weight * cmp[j].l_weight / tw, ' %'
            endif
            print, 'Weight         :', cmp[j].weight, ' +/- ', $
              cmp[j].e_weight, ' [data_unit/cmp_unit]'

            if ptr_valid(cmp[j].para) then begin
                for k = 0, n_elements(*cmp[j].para) -1 do begin
                    if strlen((*cmp[j].para)[k].name) gt 0 then begin
                        print, (*cmp[j].para)[k].name, ':', FORM='(A-15,A,$)'
                        if (*cmp[j].para)[k].dispf eq 'exp' then begin
                            val = exp((*cmp[j].para)[k].value)
                            err = val * (*cmp[j].para)[k].error
                        endif else begin
                            val = (*cmp[j].para)[k].value
                            err =  (*cmp[j].para)[k].error
                        endelse
                        if (*cmp[j].para)[k].fixed eq 0 and      $ 
                           (*cmp[j].para)[k].status eq 0 then    $
                              print, val, ' +/- ', strtrim(err, 2), $
                                     ' ', (*cmp[j].para)[k].unit $
                        else if (*cmp[j].para)[k].fixed eq 1 or  $ 
                           (*cmp[j].para)[k].status eq 1 then    $
                           print, val, ' ', (*cmp[j].para)[k].unit, ' (fixed)' $
                        else if (*cmp[j].para)[k].status eq 2 then    $
                           print, val, ' ', (*cmp[j].para)[k].unit, ' (pegged at low bound)' $
                        else if (*cmp[j].para)[k].status eq 3 then    $
                           print, val, ' ', (*cmp[j].para)[k].unit, ' (pegged at high bound)' 
                    endif
                endfor
            endif
            print,'-----------------------------------------------'

        endfor

        if nocontrib ne '' then begin
           print, 'The following component(s) had no contribution'
           print, nocontrib
           print,'-----------------------------------------------'
        endif

    endif
endfor

if readsol eq 1 then begin   ; free heap vars allocated by uly_solut_tread
    heap_free, solution.cmp
endif

end
