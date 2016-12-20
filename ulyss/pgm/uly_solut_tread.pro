;+
; NAME:
;                 ULY_SOLUT_TREAD
;
; PURPOSE:
;                 Read the solution of a fit from a .res file
;
; USAGE:
;                 solut = uly_solut_tread(<file>, STATUS=<status>)
;
; ARGUMENT:
;   <file>:       Name of the ASCII file from where to read the solution,
;                 without the '.res' extention.
;
; KEYWORD:
;   STATUS:       Returned status code, 0=OK, other value=error
;                                                            
; OUTPUT:
;   solut:        Solution structure.
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
function uly_solut_tread, file, STATUS=status

status = 0

if (size(file))[0] gt 1 or size(file,/TYPE) ne 7 then begin
   status = 1
   message, 'Invalid file name', /INFO
   print, 'Usage: solut = uly_solut_tread(<file>, STATUS=status)'
   return, 0
endif

filen = file[0]
if file_test(filen) eq 0 then filen = file + '.res'
if file_test(filen) eq 0 then begin
    status = 1
    message, 'File not found ('+file+')', /INFO
    return, 0
endif

nlines = file_lines(filen)
nread = 0  ; init counter of (solution)s read
strbuf=strarr(1)
cols_losvd = ['cz', 'e_cz', 's_cz', 'sig', 'e_sig', 's_sig', $
              'h3', 'e_h3', 's_h3', 'h4', 'e_h4', 's_h4', $
              'h5', 'e_h5', 's_h5', 'h6', 'e_h6', 's_h6']

head_sec = 0

openr, lun, filen, /GET_LUN

for j = 0L, nlines[0]-1 do begin
    readf, lun, strbuf

    if strmid(strbuf,0,7) eq '!BEGIN ' then begin ; start header section
       head_sec = 1
    endif else if head_sec gt 0 then begin ; process header recordq
       if strmid(strbuf,0,1) ne '!' then  begin
          status = 1
          message, 'Miss-formed header', /INFO
          return, 0
       endif
       if head_sec eq 1 then hdr = [strmid(strbuf,1)] $
       else hdr = [hdr, strmid(strbuf,1)]
       head_sec ++
       if strmid(strbuf,0,5) eq '!END ' then head_sec = 0

    endif else if strmid(strbuf,0,6) eq '!cols:' then begin
       if n_elements(hdr) gt 0 then begin
          sol_vers =  sxpar(hdr,'SOL_VERS', COUNT=cnt)
;         A change of the integer part of version number mark a full 
;         incompatibility (the reading pgm will not work or give wrong content)
;         A change in the decimal part indicate some incompatibility: in some
;         cases some information will not be read.
          if cnt eq 1 then begin
             if floor(sol_vers) gt 2.0 then message, /INFO, 'ULySS version not forward-compatible with this new SOLUT format, please upgrade' $
             else if sol_vers gt 2.0 then message, /INFO, 'ULySS version not fully forward-compatible with this new SOLUT format: some information may not be read, please upgrade'
          endif 
       endif
       if n_elements(sol_vers) eq 0 then sol_vers = 1

       cols = strsplit(strmid(strbuf,6),/EXTRACT)

       nc_chi = where(cols eq 'chisq', ncnt)
       if ncnt le 0 then undefine, nc_chi

       nc_snr = where(cols eq 'snr', ncnt)
       if ncnt le 0 then undefine, nc_snr

       nc_losvd = intarr(n_elements(cols_losvd))
       mask = intarr(n_elements(cols)) + 1
       for k = 0,n_elements(cols_losvd)-1 do begin
          nc_losvd[k] = where(cols eq cols_losvd[k], count)
          if count eq 1 then mask[nc_losvd[k]] = 0
       endfor

       nc_nolosvd = where(mask eq 1, count)
       if count eq 0 then message, 'Invalid solut file'

;      search the cmp columns
       nc_cmp = where(strmatch(cols[nc_nolosvd],'*_*'), count)


       if count eq 0 then message, 'Invalid solut file, no cmp information'
       nc_cmp = nc_nolosvd[nc_cmp]

;       extract the losvd columns
       ncols_losvd = where(nc_losvd ge 0, count)
       if count gt 0 then begin 
          nc_losvd = nc_losvd[ncols_losvd] 
          if sol_vers eq 1 then m = 2 else m = 3 
          moments = count / m
          if m*moments ne count then begin
             status = 1
             message, 'Invalid file format (losvd)'+strbuf, /INFO
             close, lun
             free_lun, lun
             return, 0
          endif
          losvd = dblarr(moments, /NOZERO)
          e_losvd = dblarr(moments, /NOZERO)
          s_losvd = intarr(moments, /NOZERO)
       endif else begin
          undefine, nc_losvd
          undefine, losvd
          undefine, e_losvd
          undefine, s_losvd
       endelse
        
       cols_cmp = cols[nc_cmp]
       cmpnm = gettok(cols_cmp,'_')

       idx = where(cmpnm ne shift(cmpnm,1), c)
       if c eq 0 then cmpnames = [cmpnm[0]] $
       else cmpnames = cmpnm[idx]

       n_cmp = n_elements(cmpnames)

;      cmp strucure, without the tags that we do not use here
       cmp = replicate( {name:'', $
                          para:ptr_new(), $
                          weight:0d, $
                          e_weight:0d, $
                          l_weight:0d $
                         }, n_cmp)

       if sol_vers eq 1 then m = 2 else m = 3 
       cmpn = 0
       cmpname = ''
       k = 0
       while k lt n_elements(cols_cmp) do begin
          cmp[cmpn].name = cmpnames[cmpn]
          if cols_cmp[k] ne 'weight' or cols_cmp[k+1] ne 'e_weight' or $
             cols_cmp[k+2] ne 'l_weight' then begin
             status = 1
             message, 'invalid file format (weight expected as 1st param of cmp)', /INFO
             close, lun
             free_lun, lun
             return, 0
          endif
          k += 3
          
          if k lt n_elements(cols_cmp) then begin
             while k lt n_elements(cols_cmp)-2 and cmpnm[k] eq cmpnames[cmpn] do begin
                if cols_cmp[k] eq 'weight' and cols_cmp[k+1] eq 'e_weight' and cols_cmp[k+2] eq 'l_weight' then begin
                   status = 1
                   message, 'It seems that two components have the same name ... they could not be differentiated', /INFO
                   if cmpn lt n_cmp-1 then begin
                      cmpnames = [cmpnames[0:cmpn], cmpnames[cmpn], cmpnames[cmpn+1:*]] 
                   endif else begin
                      cmpnames = [cmpnames[0:cmpn], cmpnames[cmpn]]
                   endelse
                   cmp = [cmp, {name:'', $
                                para:ptr_new(), $
                                weight:0d, $
                                e_weight:0d, $
                                l_weight:0d $
                               }]
                   n_cmp ++
                   break
                endif 
                if 'e_'+cols_cmp[k] ne cols_cmp[k+1] then begin
                   status = 1
                   message, 'Invalid file format (missing error column, col='+cols_cmp[k]+')', /INFO
                   close, lun
                   free_lun, lun
                   return, 0
                endif
                if m eq 3 then begin
                   if 's_'+cols_cmp[k] ne cols_cmp[k+2] then begin
                      status = 1
                      message, 'Invalid file format (missing status column, col='+cols_cmp[k]+')', /INFO
                      close, lun
                      free_lun, lun
                      return, 0
                   endif
                endif
                if 'g_'+cols_cmp[k] ne cols_cmp[k+m] then begin
                   status = 1
                   message, 'Invalid file format (missing guess column, col='+cols_cmp[k]+')', /INFO
                   close, lun
                   free_lun, lun
                   return, 0
                endif
                if ptr_valid(cmp[cmpn].para) eq 0 then begin
                   cmp[cmpn].para = ptr_new( $
                                    {name:cols_cmp[k], unit:'', guess:ptr_new(), step:1D-2, $
                                     limits:[0d,0d], limited:[1,1], fixed:0S, $
                                     value:0D, error:0D, status:0S, dispf:''})
                endif else begin
                   *cmp[cmpn].para = [*cmp[cmpn].para, $
                                      {name:cols_cmp[k], unit:'', guess:ptr_new(), step:1D-2, $
                                       limits:[0d,0d], limited:[1,1], fixed:0S, $
                                       value:0D, error:0D, status:0S, dispf:''}]
                endelse
                k = k + m + 1
                if k ge n_elements(cols_cmp) then break
             endwhile
          endif
          cmpn ++
       endwhile
       if n_elements(solution) eq 0 then begin
          if n_elements(losvd) gt 0 then $
             solution = {hdr:ptr_new(hdr), chisq:0., snr:0., losvd:losvd, e_losvd:e_losvd, s_losvd:s_losvd, cmp:cmp} $
          else $
             solution = {hdr:ptr_new(hdr), chisq:0., snr:0., cmp:cmp} 
       endif else begin
          if size(cmp,/DIM) ne size(solution[0].cmp,/DIM) then begin
             print, 'Record number: ', n_elements(solution)+1
             message, 'Inconsistency in this file, the number of cmp changes', /INFO
             close, lun
             free_lun, lun
             return, 0
          endif
          if n_elements(losvd) gt 0 then $
             solution = [solution, {hdr:ptr_new(hdr), chisq:0., snr:0., losvd:losvd, e_losvd:e_losvd, s_losvd:s_losvd, cmp:cmp}] $
          else $
             solution = [solution, {hdr:ptr_new(hdr), chisq:0., snr:0., cmp:cmp}]
       endelse
       
    endif else if strmid(strbuf,0,6) eq '!unts:' then begin
       if n_elements(solution) eq 0 then begin
          status = 1
          message, 'Invalid file format (no !cols before !unts)', /INFO
          return, 0
       endif
       n = n_elements(solution) - 1
;       The unts string has one item for each param + one for weight in cmps
       unts = strsplit(strmid(strbuf,6),/EXTRACT)
       
       nv = 0
       for k=0,n_elements(solution[n].cmp)-1 do begin
          nv ++                 ; yet no unit for the weight, just skip
          if ptr_valid(solution[n].cmp[k].para) then begin
             for ji=0,n_elements(*solution[n].cmp[k].para)-1 do begin
                if nv lt n_elements(unts) then begin
                   unt = strsplit(unts[nv] , ',', /EXTRACT)
                   (*solution[n].cmp[k].para)[ji].unit = unt[0] 
                   if n_elements(unt) gt 1 then $
                      (*solution[n].cmp[k].para)[ji].dispf = unt[1] 
                endif
                nv ++
             endfor
          endif
       endfor
    endif else if strmid(strbuf,0,1) ne '!' then begin
       if n_elements(solution) eq 0 then begin
          status = 1
          message, 'Invalid file format (no header)', /INFO
          close, lun
          free_lun, lun
          return, 0
       endif
       
       if n_elements(solution) eq nread then $
          solution = [solution, solution[n_elements(solution)-1]]
       
       dblbuf = dblarr(n_elements(cols), /NOZERO)
       reads, strbuf, dblbuf

       if n_elements(nc_losvd) gt 0 then begin
          if sol_vers eq 1 then m = 2 else m = 3
          losvdbuf = dblbuf(nc_losvd)
          losvdbuf = reform(losvdbuf, m, n_elements(nc_losvd)/m, /OVER)
          solution[nread].losvd = losvdbuf[0,*]
          solution[nread].e_losvd = losvdbuf[1,*]
          if m eq 3 then solution[nread].s_losvd = floor(losvdbuf[2,*])
       endif 
       
       if n_elements(nc_chi) eq 1 then solution[nread].chisq = dblbuf[nc_chi] $
       else solution[nread].chisq = 0
       
       if n_elements(nc_snr) eq 1 then solution[nread].snr = dblbuf[nc_snr] $
       else solution[nread].snr = 0
       
       if sol_vers eq 1 then m = 3 else m = 4 ; nb val for each param
       dblbuf = dblbuf[nc_cmp]
       for k=0,n_elements(solution[nread].cmp)-1 do begin
          solution[nread].cmp[k].weight = dblbuf[0]
          solution[nread].cmp[k].e_weight = dblbuf[1]
          solution[nread].cmp[k].l_weight = dblbuf[2]
          param = solution[nread].cmp[k].para
          if ptr_valid(param) then begin
             npara = n_elements(*param)
             if npara gt 0 then begin                
                pbuf = reform(dblbuf[3:2+npara*m], m, npara)
                (*param).value = reform(pbuf[0,*])
                (*param).error = reform(pbuf[1,*])
                if m eq 4 then begin
                   (*param).status = floor(reform(pbuf[2,*]))
                   (*param).fixed = ((*param).status eq [1, 1, 1, 1])
                endif
                for l=0, n_elements(*param)-1 do $
                   (*param)[l].guess = ptr_new(reform(pbuf[m-1,l]))
             endif
             next = 3 + npara * m
          endif else next = 3
          if k lt n_elements(solution[nread].cmp)-1 then dblbuf = dblbuf[next:*]
       endfor
              
       nread ++
       
    endif
 endfor
close, lun
free_lun, lun

if n_elements(solution) gt 0 then return, solution $
else return, 0

end
